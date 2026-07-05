# Where coverage lives downstream: `smaht-dac` and `parklab`

Research notes for the Park Lab take-home (per-chromosome average coverage on the
COLO829T CRAM). The point of this document is to show *who consumes a coverage
number* in the two GitHub orgs behind this role, how they compute it, and what
that implies for how the deliverable should be framed.

Everything below was traced from public code in `github.com/smaht-dac` and
`github.com/parklab` (July 2026). File paths are given so claims are checkable.

---

## TL;DR

There are **two distinct meanings of "coverage"** across these orgs, and they use
different tools for different reasons:

| | **SMaHT DAC** (`smaht-dac`) | **Park Lab** (`parklab`) |
|---|---|---|
| Tool | **mosdepth** (`-n --fast-mode`), Picard `CollectWgsMetrics` | **`samtools depth`** (default) / GATK `DepthOfCoverage` (option). **No mosdepth.** |
| Granularity | genome-wide summary + per-contig rows | **base-pair resolution**, parallelized **per chromosome** |
| Purpose | user-facing **QC metric** on the data portal | **scientific input** to variant calling + single-cell amplification QC |
| Style | production CWL/Docker → `parse-qc` → portal | Snakemake + **hand-rolled awk/R**, tabix-indexed |
| Where a human sees it | `data.smaht.org` QC card ("Estimated Average Coverage") | intermediate file feeding a caller; not surfaced as a metric |

The assignment (per-chromosome average per-base coverage) sits exactly on the
seam: it is the *mosdepth summary* that the DAC computes, produced the *samtools
depth* way that the Lab prefers.

---

## Part 1 — SMaHT DAC: coverage as a portal QC metric

### The pipeline path (tool → parser → portal)

The BAM/CRAM QC meta-workflows in **`smaht-dac/qc-pipelines`**
(`portal_objects/metaworkflows/paired-end_short_reads_BAM_quality_metrics_GRCh38.yaml`,
`samtools_mosdepth_bam.yaml`, plus long-/ultra-long-read variants) fan out to a
set of tools, each emitting a text file. A single `parse-qc` step then maps the
raw outputs to human-readable metrics and emits a QC-values JSON that becomes a
`QualityMetric` item shown on each file's **Interactive QC Assessment** page.

mosdepth is **also embedded directly in the production alignment meta-workflows**
(`smaht-dac/main-pipelines/portal_objects/metaworkflows/Illumina_alignment_GRCh38.yaml`
and the Hi-C / ONT / PacBio variants), so it runs on essentially every aligned
file, and its output feeds a `qc_ruleset` with pass/warn/fail `qc_thresholds`.

### The three coverage-producing tools and their portal metric names

Metric names, tooltips, and visibility come from
**`smaht-dac/qc-parser`** → `src/metrics_to_extract.py` (pip `qc-parser==0.8.3`):

| Tool | Invocation | Portal-facing metric |
|---|---|---|
| **mosdepth** | `mosdepth -n --fast-mode` → `.mosdepth.summary.txt` | **"Estimated Average Coverage [mosdepth]"** (from the `total` row; `visible: True`) |
| **Picard `CollectWgsMetrics`** | MQ/BQ ≥ 20, coverage cap 1000 | **"Mean Coverage (chr22) [Picard]"** + **"Coverage Standard Deviation (chr22) [Picard]"** |
| **bamstats** | — | **"Estimated Average Coverage [bamstats]"** |

```python
# smaht-dac/qc-parser : src/metrics_to_extract.py
mosdepth_metrics = {
    "total": {
        "key": "Estimated Average Coverage [mosdepth]",
        "tooltip": "Estimated average coverage",
        "derived_from": "mosdepth:total",
        "type": float,
        "visible": True,
    },
}
```

Two design tells worth noting:
- mosdepth runs with `-n --fast-mode` (skip per-base output; summary only).
- Picard's genome-wide mean is **deliberately approximated using chr22 only**
  ("approximation using chr22 only" is literally the tooltip). At WGS scale on
  AWS, fast/approximate coverage is the accepted tradeoff.

### What happens to `.mosdepth.summary.txt` itself

It is **not discarded**, but it has a narrow fate:

1. **The `total` row** → the one *visible* portal number
   ("Estimated Average Coverage [mosdepth]"). The per-chromosome rows are **not**
   parsed into visible metrics.
2. **The whole file** → the `parse-qc` step's second output is `metrics_zip`
   (`zipped: True`), fed both `MOSDEPTH_SUMMARY` and `MOSDEPTH_OUTPUT`. So the
   complete summary (all per-contig rows) plus the raw regions output are bundled
   into a QC artifact attached to the file's `QualityMetric` item — downloadable,
   but not rendered.
3. The `.mosdepth.global.dist.txt` (cumulative coverage distribution — the usual
   input to a coverage curve) is globbed and bundled too, but I found **no
   plotting / MultiQC step that visualizes it** in public code.

**No downstream *analysis* consumes it as input.** Nothing in
`smaht-dac/calling-pipelines` reads mosdepth output; the mosdepth reference in
`sentieon-pipelines/.../dnascopehybrid` is that tool's own internal coverage run,
not this file. On the DAC side, coverage is a **reporting / QC artifact
end-to-end**, not an input to variant or CNV calling.

Front-end consumers (`smaht-dac/smaht-portal`:
`types/quality_metric.py`, `static/components/viz/QualityMetricVisualizations/utils.js`,
`SubmissionStatusFileGroupQcModal.js`) read the parsed JSON, not the raw file.

---

## Part 2 — Park Lab: coverage as a variant-calling / single-cell-QC input

**Park Lab uses no mosdepth anywhere** (org-wide search is empty). Coverage code
is `samtools depth` / GATK `DepthOfCoverage`, hand-rolled, and lives inside
analysis pipelines rather than a reporting layer.

### SCAN2 — the clearest and most relevant case

**`parklab/SCAN2`** (single-cell somatic SNV calling; also `scan-snv`, `r-scan2`)
builds **base-pair-resolution depth profiles**. The header of
`snakemake/snakefile.depth_profile` is almost a direct commentary on this
assignment:

```
# Basepair resolution depth profiles are used for extrapolation of mutation burden.
# There are two tools for calculating this bp-res depth:
#       1. GATK DepthOfCoverage
#       2. samtools depths
# GATK DepthOfCoverage is ~100x slower than samtools depth and gives comparable
# output. Because of this, the two tools are parallelized differently. GATK uses
# analysis_regions and samtools depths splits jobs per chromosome.
```

The actual `samtools depth` invocation (`snakefile.depth_profile`), with the read
filters they consider correct for coverage:

```bash
samtools depth \
    -a \
    -f {input.argfile} \
    --reference {config[ref]} \
    {params.regionflag} \
    --min-MQ 60 \
    --excl-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY
```

Key choices:
- **`-a`** — report *all* positions including zero depth. Essential for a correct
  average over a whole chromosome (denominator = every base, not just covered
  bases).
- **`--min-MQ 60`** — only uniquely, confidently mapped reads contribute.
- **`--excl-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`** — drop unmapped,
  secondary/supplementary alignments, QC-fail, and PCR duplicates.
- **`--reference`** — required because inputs are CRAM.
- **per-chromosome scatter** — each chromosome is a separate job.

The output is then summarized by a **hand-rolled awk** program
(`scripts/summarize_depth_scatter.awk`), gathered per sample by
`summarize_depth_gather.R`, and stored as a bgzipped, **tabix-indexed** joint
depth matrix (`depth_profile/joint_depth_matrix.tab.gz`). The awk does not just
average — it builds a **2-D depth histogram** (single-cell depth × bulk depth),
clamped at `max_depth=500`, from which means/quantiles are derived downstream.
`r-scan2/R/depth.R` reads slices of that tabix table so memory stays bounded
across 10s–100s of cells.

Coverage here is a **scientific input**: depth profiles drive mutation-burden
extrapolation and calling sensitivity, not a dashboard.

### A second, distinct notion: binned read-start counts (amplification uniformity)

Same file has `rule binned_counts_profile`, which counts **reads that start in a
bin** (not per-base depth) using pre-computed equal-mappability bins:

```bash
samtools view -f 2 -F 3856 --min-MQ 60 {bam} {chrs} \
  | awk 'BEGIN{OFS="\t"} {print $3, $4-1, $4, NR}' \
  | bedtools intersect -sorted -c -a {bins} -b /dev/stdin -g {genome}
```

The comment explains why: read-start density "is more of an indication of
amplification uniformity and/or local copy number" and is used to compute MAPD.
This is a good reminder that "coverage" splinters into several definitions once
you care about *why* you're measuring it.

### Other Park Lab consumers (all samtools-based, same flavor)

- **PaSDqc** — single-cell whole-genome-amplification QC; models coverage
  *uniformity* (spectral density) to flag bad amplification.
- **MosaicForecast** — `samtools depth` for read-level features in mosaic calling.
- **LiRA**, **NGSCheckMate** (`samtools depth` at SNP sites for sample identity),
  **HiScanner**, **luquette-glia-analysis** (`digest_chrom_depth.R`) — each
  hand-rolls depth summarization for its own analysis.

---

## Part 2.5 — Tool comparison: what to use when

`chromcov` (this repo) is a finite-difference / event-accumulation depth calc over
`read.get_blocks()` (M/=/X ref-matched bases only), per chromosome, with a
samtools-style bitmask `ReadFilter` and two modes (per-base depth array vs an O(1)
aggregate). Algorithmically it is a focused reimplementation of mosdepth's core.

| Tool | What it computes | Method / what it counts | Overlap-aware? | Speed | Why you'd reach for it |
|---|---|---|---|---|---|
| **`chromcov` (this repo)** | per-chrom mean + optional per-base depth array + breadth-at-depth | pysam; **finite-difference** (+1/−1 at `get_blocks()` edges → cumsum); counts M/=/X ref-matched bases only | **No** (double-counts mate overlaps — flagged in the code) | Fast; aggregate mode is O(1) memory | You want *exactly* per-chrom average + a depth distribution, in pure Python, with fully auditable logic and no external binary. Best when the deliverable *is* the coverage logic. |
| **mosdepth** | per-base or per-region/windowed depth, per-chrom summary, cumulative dist | Same finite-difference algorithm, in Nim; `.summary.txt` (per-contig mean) + `.global.dist.txt` | **Yes, by default** (and CIGAR-aware). Both are **disabled by `--fast-mode`** — which is exactly how the DAC runs it (`mosdepth -n --fast-mode`) | Very fast (compiled, threaded); faster still in fast-mode | Production WGS at scale. The standard. What the SMaHT DAC runs. Reach for it when you need speed + a battle-tested standard and don't need to own the code. |
| **`samtools coverage`** | per-chrom table: meandepth, breadth (% covered), meanbaseq, meanmapq | pileup engine, one row per reference | Yes (pileup-level) | Fast | The fastest *one-liner* answer to this exact assignment. Reach for it to sanity-check numbers or when a shell command suffices. |
| **`samtools depth`** | per-base depth (optionally all positions with `-a`) | pileup engine, streams position → depth | Configurable (base-level dedup possible) | Moderate (per-base output is heavy) | When you need raw per-base depth to pipe into custom summarization (this is what SCAN2 uses). |
| **`samtools idxstats`** | mapped **read counts** per contig | reads the `.crai/.bai` index only — **no depth scan** | N/A | Instant | Crude coverage *estimate* (`reads×readlen/contiglen`) with zero I/O over alignments. Use when "roughly how deep, per chromosome?" is enough. Ignores clipping, indels, MAPQ. |
| **`samtools` / pysam `pileup`** | per-position column of overlapping reads | builds the actual pileup column; full base/qual context | Yes (you see each fragment) | Slow | When you need *base-level* detail (allele counts, base quality) — i.e., variant calling, not just depth. Overkill for mean coverage. |
| **Park Lab — SCAN2 `depth_profile`** | **joint 2-D depth histogram** (single-cell depth × matched-bulk depth), per chrom, clamped at 500 | `samtools depth -a --min-MQ 60 --excl-flags …` → **hand-rolled awk** → tabix-indexed joint matrix | inherits `samtools depth` | Moderate; per-chrom parallel | Feeds SCAN2's sensitivity / burden model. See "why hand-roll" below. |
| **Park Lab — SCAN2 `binned_counts`** | **read-start density** in equal-mappability bins → MAPD | `samtools view -f2 -F3856 -q60` → bedtools intersect into pre-computed bins | N/A (counts read starts, not depth) | Moderate | Amplification uniformity + local copy number for single cells — a bespoke metric, not "coverage." |
| **Park Lab — PaSDqc** | spectral density of coverage uniformity | FFT / periodogram over the coverage signal | N/A | Moderate | QC of whole-genome-amplification artifacts in single cells. Again the *shape* of coverage, not mean depth. |

### Why the Park Lab hand-rolls its own

They do **not** hand-roll *mean coverage* — that's trivial and they'd just use
samtools. They hand-roll because the thing they actually need is not a coverage
number at all; it is a bespoke, distribution-shaped statistic no general tool
emits:

1. **The output is a joint distribution, not a scalar.** SCAN2 needs, per base,
   the depth in the single cell **and** its matched bulk simultaneously — a 2-D
   histogram (SC depth × bulk depth). mosdepth / `samtools coverage` give one
   per-sample mean; they can't produce the *paired* distribution that feeds
   SCAN2's allele-balance and sensitivity model. The awk is the novel part;
   `samtools depth` is just the feedstock.
2. **It's a model input, not a QC report.** DAC coverage is a number a human reads
   on a portal. SCAN2 coverage is a covariate inside a statistical caller, driving
   mutation-burden extrapolation and per-region sensitivity. That needs the full
   clamped distribution, tabix-indexed and sliceable per region, memory-bounded
   across 10s–100s of cells (`r-scan2/R/depth.R`). No off-the-shelf tool packages
   depth that way.
3. **Read-filter consistency with the caller.** The depth model must count reads
   *identically* to the SNV caller (MQ60, exact flag exclusions). Owning the code
   guarantees the depth profile and the variant model never drift apart.
4. **They compose, they don't reinvent.** They *do* use `samtools depth` and
   `bedtools`; the hand-rolling is confined to the glue (awk histogram, R digest).
   They also explicitly benchmarked `samtools depth` vs GATK DepthOfCoverage
   (comparable output, ~100× faster) before choosing — deliberate tool selection
   plus a thin custom layer, not NIH syndrome.
5. **`binned_counts` is the cleanest example.** It measures read-start density in
   mappability-normalized bins to compute MAPD (amplification uniformity / local
   CN) — a single-cell-specific quantity with no general-purpose tool, so it's
   hand-built.

**One-line framing:** the DAC treats coverage as a *metric* and uses the standard
tool (mosdepth); the Lab treats coverage as a bespoke *model input* and builds
exactly the distribution its callers need — reusing samtools for the heavy lifting
and hand-rolling only the novel statistic. `chromcov` sits between them:
mosdepth's algorithm reimplemented for full transparency, which is the right call
when the coverage logic itself is the deliverable.

---

## Part 3 — Implications for the deliverable

1. **Your instinct was right**: the Lab both hand-rolls *and* uses samtools, in
   different repos for different reasons. Coverage is not one thing here.

2. **The samtools-based, per-chromosome, fast framing is exactly how this group
   thinks.** SCAN2's explicit "samtools is ~100× faster than GATK, so we
   parallelize per chromosome" is independent confirmation.

3. **Read filters matter and are worth stating.** SCAN2's
   `--min-MQ 60 --excl-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY` and `-a`
   (count zero-depth positions) are the defensible defaults. Whatever your tool
   does, document the filter and the denominator explicitly — that is the part a
   reviewer will scrutinize.

4. **Know both audiences.** To the DAC, "coverage" is a *QC metric* surfaced on a
   portal with a stable key + tooltip + machine-readable JSON. To the Lab, it is a
   *calling input* / burden-extrapolation profile. A one-line note that you
   understand both sides shows range.

5. **COLO829T has CNVs**, so the *per-chromosome* breakdown is the interesting
   part — chromosome/arm means deviate from the genome mean in proportion to copy
   number. That's a natural sanity check to call out, and it's the exact detail
   the DAC portal keeps (in the zip) but doesn't display.

---

## Sources

- `github.com/smaht-dac/qc-pipelines` — `portal_objects/metaworkflows/*`, `descriptions/{mosdepth,parse-qc_BAM_Samtools_Mosdepth,picard_CollectWgsMetrics}.cwl`, `dockerfiles/parseqc/Dockerfile`
- `github.com/smaht-dac/qc-parser` — `src/metrics_to_extract.py`
- `github.com/smaht-dac/main-pipelines` — `portal_objects/metaworkflows/Illumina_alignment_GRCh38.yaml`
- `github.com/smaht-dac/smaht-portal` — `types/quality_metric.py`, `.../QualityMetricVisualizations/utils.js`
- `github.com/parklab/SCAN2` — `snakemake/snakefile.depth_profile`, `scripts/summarize_depth_scatter.awk`
- `github.com/parklab/r-scan2` — `R/depth.R`
- `github.com/parklab/{PaSDqc,MosaicForecast,NGSCheckMate,LiRA,HiScanner,luquette-glia-analysis}`
- SMaHT Data Portal — `https://data.smaht.org/`
- SMaHT pipelines docs — `https://smaht-dac.github.io/pipelines-docs/`
