# Coverage in Python: reference snippets

Companion to [`coverage-downstream-research.md`](./coverage-downstream-research.md).

These are **illustrative reference implementations** of the coverage patterns the
DAC and the Park Lab use, written in Python so they're easy to compare against a
hand-rolled `chromcov`. They are not copied from the repos; they reproduce the
*semantics* those repos rely on (read filters, the `-a` denominator, mosdepth's
summary format). Each is self-contained and uses `pysam` unless noted.

> **Note:** treat these as pseudo-code, not runnable code. In particular snippets
> 2 and 5 pass a `flag_filter` kwarg to `pysam.pileup` that does not exist in the
> real API — they illustrate the intended filter semantics only. Snippet 3's
> `read_callback` form is the runnable way to filter reads.

> Definitions matter. "Average per-base coverage of a chromosome" =
> `sum(depth at every base) / chromosome_length`. The subtleties are (a) which
> reads count, and (b) whether zero-depth bases are in the denominator. The
> snippets make both explicit.

---

## 1. Parse `mosdepth`'s `.mosdepth.summary.txt` (the DAC downstream format)

This is the exact file the SMaHT DAC produces and archives. Parsing it lets you
report the same per-chromosome table the DAC keeps in its QC zip, and reproduce
the single number it surfaces ("Estimated Average Coverage [mosdepth]" = the
`total` row's `mean`).

```python
from pathlib import Path

def parse_mosdepth_summary(path: str) -> dict[str, float]:
    """Return {contig: mean_depth} from a *.mosdepth.summary.txt file.

    Columns: chrom, length, bases, mean, min, max.
    mosdepth emits per-contig rows, matching `_region` rows when run with --by,
    and a final `total` row = the genome-wide mean the DAC portal displays.
    """
    means: dict[str, float] = {}
    for i, line in enumerate(Path(path).read_text().splitlines()):
        if i == 0:  # header
            continue
        chrom, length, bases, mean, mn, mx = line.split("\t")
        means[chrom] = float(mean)
    return means


summary = parse_mosdepth_summary("COLO829T_TEST.mosdepth.summary.txt")
print("Genome-wide (portal 'Estimated Average Coverage'):", summary["total"])
for chrom in [c for c in summary if c not in ("total",) and "_region" not in c]:
    print(f"{chrom}\t{summary[chrom]:.2f}")
```

---

## 2. Per-chromosome mean depth the SCAN2 way (pysam `pileup`)

Reproduces `samtools depth -a --min-MQ 60 --excl-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`
semantics: count every base (`-a`), only high-MQ reads, drop dup/secondary/etc.
This is the definition the Park Lab treats as correct for coverage. Streaming, so
memory is O(1) per chromosome.

```python
import pysam

# samtools --excl-flags UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY = 0x704
EXCLUDE = (
    pysam.FUNMAP        # 0x4
    | pysam.FSECONDARY  # 0x100
    | pysam.FQCFAIL     # 0x200
    | pysam.FDUP        # 0x400
    | pysam.FSUPPLEMENTARY  # 0x800
)
MIN_MQ = 60

def mean_depth_per_chrom(cram: str, reference: str) -> dict[str, float]:
    """Average per-base depth per chromosome, SCAN2-style filters.

    Denominator is the FULL contig length (zero-depth bases included), matching
    `samtools depth -a`.
    """
    af = pysam.AlignmentFile(cram, "rc", reference_filename=reference)
    out: dict[str, float] = {}
    for contig, length in zip(af.references, af.lengths):
        total_depth = 0
        for col in af.pileup(
            contig,
            stepper="nofilter",       # we apply our own flag/MQ filter below
            min_mapping_quality=MIN_MQ,
            flag_filter=EXCLUDE,
            truncate=True,
        ):
            total_depth += col.nsegments
        out[contig] = total_depth / length if length else 0.0
    return out
```

> `pileup` is convenient but slow on WGS. For a real tool prefer `count_coverage`
> (snippet 3) or shell out to `samtools depth`. The point here is to show the
> *filter semantics* line up with SCAN2.

---

## 3. Faster: `count_coverage` (array-based, still `-a` semantics)

`AlignmentFile.count_coverage` returns per-base A/C/G/T counts for a region using
a fast C loop. Summing the four channels gives depth; dividing by contig length
gives the same average as snippet 2, but much faster.

```python
import pysam

def mean_depth_per_chrom_fast(cram: str, reference: str) -> dict[str, float]:
    af = pysam.AlignmentFile(cram, "rc", reference_filename=reference)
    out: dict[str, float] = {}
    for contig, length in zip(af.references, af.lengths):
        # quality_threshold=0 keeps base-quality out of it; we only filter reads.
        a, c, g, t = af.count_coverage(
            contig,
            read_callback=lambda r: (
                r.mapping_quality >= 60
                and not (r.is_unmapped or r.is_secondary or r.is_qcfail
                         or r.is_duplicate or r.is_supplementary)
            ),
            quality_threshold=0,
        )
        total = sum(a) + sum(c) + sum(g) + sum(t)
        out[contig] = total / length if length else 0.0
    return out
```

---

## 4. The "fast approximation" the DAC leans on (idxstats-style, no depth scan)

Picard's chr22-only shortcut and mosdepth's `--fast-mode` both trade exactness
for speed. The cheapest estimate needs no per-base scan at all: it divides total
aligned bases on a contig by the contig length. `samtools idxstats` gives mapped
*read counts* per contig; multiply by mean read length for a bases estimate — or
sum actual aligned lengths in one pass for a better one.

```python
import pysam

def approx_coverage_per_chrom(cram: str, reference: str) -> dict[str, float]:
    """~coverage = aligned bases on contig / contig length. One linear pass,
    no per-base accumulation. Close to mosdepth --fast-mode in spirit."""
    af = pysam.AlignmentFile(cram, "rc", reference_filename=reference)
    lengths = dict(zip(af.references, af.lengths))
    aligned_bases: dict[str, int] = {c: 0 for c in af.references}
    for read in af.fetch(until_eof=True):
        if (read.is_unmapped or read.is_secondary or read.is_qcfail
                or read.is_duplicate or read.is_supplementary
                or read.mapping_quality < 60):
            continue
        aligned_bases[read.reference_name] += read.query_alignment_length
    return {c: aligned_bases[c] / lengths[c] if lengths[c] else 0.0
            for c in af.references}
```

> Caveat to note in a writeup: this counts each read's aligned length once and
> ignores CIGAR deletions/insertions distribution, so it's an *estimate*. Exact
> per-base depth (snippets 2/3) and this estimate diverge in indel-rich or
> soft-clipped regions — a good thing to mention when you report both.

---

## 5. Depth *distribution*, not just the mean (the awk histogram idea)

SCAN2's `summarize_depth_scatter.awk` builds a clamped depth **histogram**, not a
scalar mean — because for a CNV cell line like COLO829T the *shape* of the depth
distribution per chromosome is more informative than the average (a chromosome
with a copy-number gain shifts the whole distribution). This reproduces the
single-sample version and derives mean + median from it.

```python
import numpy as np
import pysam

def depth_histogram_per_chrom(cram: str, reference: str, max_depth: int = 500):
    """Per-chromosome clamped depth histogram (SCAN2-style), plus mean/median.

    Returns {contig: {"hist": np.ndarray[max_depth+1], "mean": float,
                       "median": float}}. Depth is clamped at max_depth like the
    awk (values >= max_depth land in the top bin)."""
    af = pysam.AlignmentFile(cram, "rc", reference_filename=reference)
    results = {}
    for contig, length in zip(af.references, af.lengths):
        hist = np.zeros(max_depth + 1, dtype=np.int64)
        covered = 0
        for col in af.pileup(contig, stepper="nofilter",
                             min_mapping_quality=60, truncate=True,
                             flag_filter=0x704):
            d = min(col.nsegments, max_depth)
            hist[d] += 1
            covered += 1
        # zero-depth bases (the `-a` positions pileup never emits):
        hist[0] += length - covered
        depths = np.arange(max_depth + 1)
        total_bases = hist.sum()
        mean = float((depths * hist).sum() / total_bases) if total_bases else 0.0
        cum = np.cumsum(hist)
        median = float(np.searchsorted(cum, total_bases / 2))
        results[contig] = {"hist": hist, "mean": mean, "median": median}
    return results
```

---

## Source provenance and citations

The snippets reproduce the *semantics* of these upstream files. "Last modified" is
the date of the most recent commit touching each file (default branch), retrieved
via the GitHub API in July 2026.

| Snippet(s) | Reproduces | Source file | Last modified |
|---|---|---|---|
| 1 | mosdepth summary format + DAC metric | `smaht-dac/qc-pipelines` `descriptions/mosdepth.cwl` | 2024-10-17 |
| 1 | portal metric key ("Estimated Average Coverage") | `smaht-dac/qc-parser` `src/metrics_to_extract.py` | 2025-08-07 |
| 2–5 | `samtools depth -a` filters, per-chrom scatter | `parklab/SCAN2` `snakemake/snakefile.depth_profile` | 2025-07-27 |
| 5 | clamped depth histogram | `parklab/SCAN2` `scripts/summarize_depth_scatter.awk` | 2024-06-15 |
| 2, 5 | tabix-sliced per-sample depth reads | `parklab/r-scan2` `R/depth.R` | 2025-08-22 |

### Papers that use / describe this code

- **SCAN2** (snippets 2–5) — Luquette, L.J., Miller, M.B., Kim, Z. *et al.*
  "Single-cell genome sequencing of human neurons identifies somatic point
  mutation and indel enrichment in regulatory elements." *Nature Genetics* **54**,
  1564–1571 (2022). DOI: **10.1038/s41588-022-01180-2**
  - Predecessor method / protocol (SCAN-SNV): "Somatic Single-Nucleotide Variant
    Calling from Single-Cell DNA Sequencing Data Using SCAN-SNV," *Methods in
    Molecular Biology* (2022). DOI: **10.1007/978-1-0716-2293-3_17**
- **mosdepth** (snippet 1) — Pedersen, B.S. & Quinlan, A.R. "Mosdepth: quick
  coverage calculation for genomes and exomes." *Bioinformatics* **34**(5),
  867–868 (2018). DOI: **10.1093/bioinformatics/btx699**
  *(Quinlan lab tool; used as a dependency by the SMaHT DAC QC pipelines, not
  authored by the Park Lab.)*

---

## Choosing among these for the assignment

- **Correctness reference / "the Park Lab definition":** snippet 2 or 3
  (`-a`, MQ≥60, standard flag exclusions). Report the filters explicitly.
- **Speed at WGS scale:** snippet 3 (`count_coverage`) or shell out to
  `samtools depth` and stream; snippet 4 if an approximation is acceptable.
- **Most informative report for a CNV sample:** snippet 5 — per-chromosome
  mean **and** median (or a small distribution), which surfaces the copy-number
  structure the mean alone hides.
- **Interop with the DAC:** snippet 1, so your numbers can be checked against a
  mosdepth run and reported in the same shape the portal archives.
