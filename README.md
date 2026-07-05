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
Full combined table (mean + stats + copy number + flags), all chromosomes:
`out/<coverage-key>/<analysis-key>/coverage.tsv`.

## Quick start

```bash
# 1. get the data (see dev/CONTEXT.md for URLs) into data/
#    COLO829T_TEST.cram(.crai) + GCA_000001405.15_GRCh38_no_alt_analysis_set.fa(.fai)
uv sync                         # installs the package + the `chromcov` CLI

# full run (default): per-base depth -> combined stats table + windows + copy
# number + strata + plots + reusable tracks, all written under ./out
chromcov coverage --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# just the numbers: mean-only per-chromosome table (no per-base, stats, or plots)
chromcov coverage --cram … --reference … --fast

# whole genome by default; add `--chroms chr21,chrM` for a fast subset,
# `--config config.example.yml` to drive a run from one file, or
# `--strata easy=SMaHT_easy_hg38.bed.gz,...` for callability-stratified coverage.
```

The combined coverage table (mean, median, MAD, breadth, copy number, flags) is
one file — `--fast` just leaves the per-base columns unfilled. `--write` also
archives the table under `runs/<name>/`; `chromcov collate` compares them. An
optional mosdepth cross-check (`scripts/mosdepth_coverage.py`) emits the same
`coverage.tsv` format for validation — it needs the `mosdepth` binary on PATH,
but chromcov itself does not.

## Compute once, re-analyze cheaply (per-base tracks + reuse)

Coverage is expensive; re-deriving stats/strata/plots from it is cheap. So the
per-base depth is a first-class, interoperable **output** (a standard
`bedgraph.gz`, one per chromosome), written by default — not a hidden cache; see
`WORKFLOW.md`. Three levels, each content-addressed:

```bash
chromcov fetch strata                    # download the SMaHT easy/difficult/extreme BEDs

# 1st run writes per-base tracks (Level 1) under out/<coverage-key>/
chromcov coverage --cram … --reference … --chroms chr20,chr21,chrX,chrY,chrM

# a 2nd run REUSES those tracks (no CRAM recompute) and nests a separate hashed
# run dir beside them, so stratified-vs-not compares off one coverage dataset
chromcov coverage --cram … --reference … --chroms chr20,chr21,chrX,chrY,chrM \
  --strata easy=data/SMaHT_easy_hg38.bed.gz,difficult=data/SMaHT_difficult_hg38.bed.gz,extreme=data/SMaHT_extreme_hg38.bed.gz
```

The Level-1 key hashes *inputs + read-filter params* only, so changing `--window`
or `--strata` reuses the tracks and re-runs just the cheap reductions. Everything
for one coverage dataset is co-located:

```
out/<coverage-key>/                       chrN.per-base.bedgraph.gz + coverage.json   (Level 1)
out/<coverage-key>/<analysis-key>/        coverage.tsv · windows · strata · plots + run.json  (Level 2)
```

The same steps run under **Snakemake** (scatter one track per chromosome, gather
into one full run; free parallelism + resume):

```bash
uv run snakemake --cores 4 --configfile config/config.example.yaml        # -n for a dry-run DAG
```

## What it computes

**Core (the deliverable):** per-chromosome mean coverage, via an event-based
finite-difference algorithm that costs O(reads), not O(bases).

**Extensions (what a production QC pipeline wants):**
- **Robust stats** — median, scaled MAD, robust CV (MAD/median), IQR, breadth-at-depth. Robust measures matter because multi-mapping pileups make plain sd/mean useless (chr21: sd 71 vs MAD 7).
- **Optional mosdepth cross-check** — a standalone add-on (`scripts/mosdepth_coverage.py`) runs `mosdepth` and converts its output to the same `coverage.tsv` format, so validation is an explicit `diff` rather than a second backend in the core.
- **Callability stratification** — coverage within Park Lab's SMaHT easy/difficult/extreme region tiers; the "easy" coverage is the variant-callable number worth reporting.
- **Approximate copy number** — depth normalized to the callability-masked autosomal median; windowed CN scatter reveals intrachromosomal breakpoints.
- **Abnormality flags** — aneuploidy-aware CN gain/loss/depletion, low-depth, uneven, low-callability, extreme-depth.
- **Reproducibility** — input preflight (sorted / indexed / reference-M5 match), provenance sidecars, Docker + CWL sketches.

## Code flow

```
chromcov coverage ──► preflight (validate.py)   sorted? indexed? reference matches?
       │
       ├── --fast ─► dispatch.run_coverage ─► calc_cov ─► mean-only coverage.tsv
       │             (optional cross-check: scripts/mosdepth_coverage.py, same format)
       │
       └── default ─► CoverageAnalysis (per-base, one memory-bounded pass per chrom):
                 calc_cov(per_base=True) ──► ChromDepth(per-base depth vector)
                     ├─ .histogram() → DepthHistogram ─► stats / MAD / breadth
                     ├─ .windowed_means() ─► windows ─► plots, copy number
                     ├─ Strata.mask()     ─► per-tier callable coverage
                     ├─ qc.chrom_flags    ─► abnormality flags → combined coverage.tsv
                     └─ PerBaseStore      ─► per-base bedgraph tracks (reused next run)
```

## Repo layout

| Path | What |
|---|---|
| `chromcov/` | the installable package (config, coverage, analysis, per-base tracks, CLI) |
| `scripts/mosdepth_coverage.py` | optional mosdepth cross-check add-on (same output format) |
| `tests/` | unit tests (reduction/QC math, track round-trip) + native↔mosdepth cross-check |
| `WORKFLOW.md` | mermaid diagrams of the three-level workflow + orchestration layers |
| `workflow/Snakefile`, `config/` | Snakemake orchestration layer + example config |
| `config.example.yml` | annotated run-config driving `chromcov --config` |
| `dev/reproducibility-sketch/` | CWL tool + Dockerfile (call `chromcov coverage`) |
| `dev/CONTEXT.md` | the assignment brief |
| `DEVELOPMENT.md` | algorithm + design deep dive |
| `TODO.md` | what remains to make this a finished, sendable repo |
