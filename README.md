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
uv sync

# 2. per-chromosome coverage table (native backend, no container needed)
uv run python dev/dispatch-sketch/run_example.py --backend native --write

# 3. full analysis: stats + windows + copy number + plots + QC flags
uv run --with matplotlib python dev/dispatch-sketch/analyze.py
```

Optional callability stratification (Park Lab SMaHT masks) and the mosdepth
cross-check are described in `dev/dispatch-sketch/README.md`.

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
       ├──► dispatch.run_coverage ──► native.py  ─┐  both emit
       │                          └► mosdepth.py ─┴► list[ChromCoverage]
       │                                              └► output.py (run dir + provenance)
       │
       └──► analyze.py (native per-base, one memory-bounded pass per chrom):
                 calc_cov(per_base=True) ──► per-base depth vector
                     ├─ analysis.depth_histogram ─► stats / MAD / breadth
                     ├─ analysis.windowed_means  ─► windows ─► plots, copy number
                     ├─ strata.stratum_mask      ─► per-tier callable coverage
                     ├─ qc.chrom_flags           ─► abnormality flags
                     └─ analysis.rle_intervals   ─► optional per-base bedgraph
```

## Repo layout

| Path | What |
|---|---|
| `chromcov/` | the installable package (**in-progress** — see `TODO.md`) |
| `dev/dispatch-sketch/` | the working reference implementation (runnable sketches) |
| `dev/reproducibility-sketch/` | provenance, CWL tool, Dockerfile |
| `dev/CONTEXT.md` | the assignment brief |
| `DEVELOPMENT.md` | algorithm + design deep dive |
| `TODO.md` | what remains to make this a finished, sendable repo |

> **Current state, honestly:** the polished code lives as runnable *sketches*
> under `dev/dispatch-sketch/`. Consolidating them into the installable
> `chromcov` package with a CLI, real tests, and CI is the top of `TODO.md`.
