# chromcov — per-chromosome coverage from a CRAM

Computes **average per-base sequence coverage for each chromosome** in a CRAM,
and reports it in the ways a bioinformatics analyst and a software engineer would
each want: a clean per-chromosome table, robust dispersion stats, callability-
stratified coverage, an approximate copy-number readout, plots, and full run
provenance.

> **The task** (Park Lab / SMaHT DAC take-home): *"Develop a program that
> calculates the average per-base sequence coverage for each chromosome in the
> CRAM file we are sharing."* Sample: **COLO829T**, a CNV-rich benchmark cancer
> line, aligned to **GRCh38** (`no_alt_analysis_set`). See `dev/CONTEXT.md`.

## The answer

Average per-base coverage = **aligned bases / chromosome length**. For the
provided COLO829T CRAM (~15× WGS), the primary assembly:

| chrom | mean | flag | chrom | mean | flag |
|---|---:|---|---|---:|---|
| chr1 | 15.2 | | chr13 | 13.8 | |
| chr3 | 19.6 | gain | chr20 | 22.5 | gain |
| chr5 | 10.5 | loss | chrX | 10.0 | (1 copy, male) |
| chr7 | 21.2 | gain | chrY | 0.8 | **loss of Y** |
| chr10 | 11.0 | loss | chrM | 19790 | (mitochondrial) |

The per-chromosome spread is real COLO829T biology — aneuploidy (chr3/7/20 gains,
chr5/10/18 losses), a single X, loss of Y, and high mitochondrial copy number.
Full table + all chromosomes: `out/coverage.stats.tsv`.

## Quick start

```bash
# 1. get the data (see dev/CONTEXT.md for URLs) into data/
#    COLO829T_TEST.cram(.crai) + GCA_000001405.15_GRCh38_no_alt_analysis_set.fa(.fai)
uv sync                         # installs the package + the `chromcov` CLI

# 2. per-chromosome coverage table (the deliverable; native backend, no container)
chromcov coverage --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# 3. full analysis (superset of the table): stats + windows + copy number + plots + QC flags
chromcov analyze --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# whole genome by default; add `--chroms chr21,chrM` for a fast subset,
# `--config config.example.yml` to drive a run from one file, or
# `--strata easy=SMaHT_easy_hg38.bed.gz,...` for callability-stratified coverage.
```

`chromcov coverage --write` archives each run under `runs/<name>/` with a
provenance sidecar; `chromcov collate` compares them. The mosdepth cross-check
backend (`--backend mosdepth`) needs the binary on PATH.

## Compute once, analyze many (per-base tracks + reuse)

Coverage is expensive; re-deriving stats/strata/plots from it is cheap. So the
per-base depth is a first-class, interoperable **output** (a standard
`bedgraph.gz`, one per chromosome), not a hidden cache — see `WORKFLOW.md` for
the diagram. Three levels, each content-addressed:

```bash
chromcov fetch strata                    # download the SMaHT easy/difficult/extreme BEDs

# 1st analyze --per-base writes per-base tracks (Level 1) under out/perbase/<coverage-key>/
chromcov analyze --cram … --reference … --chroms chr20,chr21,chrX,chrY,chrM --per-base

# a 2nd analyze REUSES those tracks (no CRAM recompute) and writes a separate
# hashed run dir, so stratified-vs-not is a fast comparison off the same coverage
chromcov analyze --cram … --reference … --chroms chr20,chr21,chrX,chrY,chrM \
  --strata easy=data/SMaHT_easy_hg38.bed.gz,difficult=data/SMaHT_difficult_hg38.bed.gz,extreme=data/SMaHT_extreme_hg38.bed.gz
```

The Level-1 key hashes *inputs + read-filter params* only, so changing `--window`
or `--strata` reuses the tracks and re-runs just the cheap reductions. Each
analysis lands in `out/analysis/<analysis-key>/` with a `run.json` sidecar.

The same steps run under **Snakemake** (scatter one track per chromosome, gather
into one analysis; free parallelism + resume):

```bash
uv run snakemake --cores 4 --configfile config/config.example.yaml        # -n for a dry-run DAG
```

## What it computes

**Core (the deliverable):** per-chromosome mean coverage, via an event-based
finite-difference algorithm that costs O(reads), not O(bases).

**Extensions (what a production QC pipeline wants):**
- **Robust stats** — median, scaled MAD, robust CV (MAD/median), IQR, breadth-at-depth. Robust measures matter because multi-mapping pileups make plain sd/mean useless (chr21: sd 71 vs MAD 7).
- **Two interchangeable backends** — a hand-rolled `pysam` calculator and `mosdepth`, behind one config, so each validates the other.
- **Callability stratification** — coverage within Park Lab's SMaHT easy/difficult/extreme region tiers; the "easy" coverage is the variant-callable number worth reporting.
- **Approximate copy number** — depth normalized to the callability-masked autosomal median; windowed CN scatter reveals intrachromosomal breakpoints.
- **Abnormality flags** — aneuploidy-aware CN gain/loss/depletion, low-depth, uneven, low-callability, extreme-depth.
- **Reproducibility** — input preflight (sorted / indexed / reference-M5 match), provenance sidecars, Docker + CWL sketches.

## Code flow

```
CoverageConfig ──► preflight (validate.py)   sorted? indexed? reference matches?
       │
       ├──► dispatch.run_coverage ──► NativeBackend  ─┐  both emit
       │                          └► MosdepthBackend ─┴► list[ChromCoverage]
       │                                                  └► RunStore (run dir + provenance)
       │
       └──► CoverageAnalysis (native per-base, one memory-bounded pass per chrom):
                 calc_cov(per_base=True) ──► ChromDepth(per-base depth vector)
                     ├─ .histogram()      ─► stats / MAD / breadth
                     ├─ .windowed_means() ─► windows ─► plots, copy number
                     ├─ Strata.mask()     ─► per-tier callable coverage
                     ├─ qc.chrom_flags    ─► abnormality flags
                     └─ .rle_intervals()  ─► optional per-base bedgraph
```

## Repo layout

| Path | What |
|---|---|
| `chromcov/` | the installable package (config, backends, analysis, per-base tracks, CLI) |
| `tests/` | unit tests (reduction/QC math, track round-trip) + native↔mosdepth cross-check |
| `WORKFLOW.md` | mermaid diagrams of the three-level workflow + orchestration layers |
| `workflow/Snakefile`, `config/` | Snakemake orchestration layer + example config |
| `config.example.yml` | annotated run-config driving `chromcov --config` |
| `dev/reproducibility-sketch/` | CWL tool + Dockerfile (call `chromcov coverage`) |
| `dev/CONTEXT.md` | the assignment brief |
| `DEVELOPMENT.md` | algorithm + design deep dive |
| `TODO.md` | what remains to make this a finished, sendable repo |
