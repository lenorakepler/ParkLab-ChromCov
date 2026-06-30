"""Generate a comparison of the coverage tools explored in this project,
as both Markdown (docs/coverage_tools_comparison.md) and a styled HTML page
(docs/coverage_tools_comparison.html). Single source of truth below.
"""
import re
import html
from pathlib import Path

DOCS = Path().resolve() / "docs"
DOCS.mkdir(exist_ok=True)

COLUMNS = ["Tool (engine)", "Run via (this project)", "CRAM ref", "Speed (17 GB CRAM)",
           "Main outputs", "Binning / granularity", "Key caveats", "Best use case"]

# Practical, interchangeable coverage tools (all run this session via biocontainers)
ROWS = [
	["mosdepth (C / htslib)",
	 "**Docker biocontainer** `quay.io/biocontainers/mosdepth:0.3.3--h37c5b7d_2`<br>`mosdepth --fasta REF --by N PREFIX CRAM`<br>_(also: bioconda, precompiled binary)_",
	 "`--fasta`",
	 "Fast (streams; minutes)",
	 "`.summary.txt` per-chrom mean; `.regions.bed.gz` windows; `.global/region.dist.txt` breadth-at-depth; `.thresholds.bed.gz`; d4",
	 "Native `--by` windows/BED; `--thresholds`; `--quantize`",
	 "Needs ref for CRAM; mate-overlap aware (≈7% lower abs. depth than samtools/pandepth); per-base mode slower",
	 "The standard: windowed depth, QC breadth distributions, CNV binning; production / clinical pipelines"],

	["PanDepth (C++)",
	 "**Docker biocontainer** `quay.io/biocontainers/pandepth:2.26--h02cc909_0`<br>`pandepth -i CRAM -r REF [-w N|-b BED|-g GFF] -o PREFIX`<br>_(also: bioconda, GitHub Linux binary, source)_",
	 "`-r`",
	 "Fastest streaming tool",
	 "`.chr.stat.gz` per-chrom mean + cov%; `.win.stat.gz` windows; GC with `-c`",
	 "`-w` windows, `-b` BED, `-g` GFF",
	 "No mate-overlap correction (≈7% high vs mosdepth); breadth at one `-d` threshold/run; newer, less validated",
	 "Fast per-chrom + windowed coverage, GC content; quick CNV exploration; reads CRAM directly"],

	["samtools coverage (htslib / pysam)",
	 "**Python / pysam** in the uv venv — no Docker (pysam wheel bundles htslib)<br>`pysam.coverage('--reference', REF, '-r', chrom, CRAM)`<br>_(or samtools CLI)_",
	 "`--reference`",
	 "Streams; slow when looped per-chrom",
	 "Table: numreads, covbases, coverage%, meandepth, meanbaseq, meanmapq",
	 "Per-region (`-r`); `-m` ASCII histogram; no windowing",
	 "Breadth only at ≥1×; counts overlapping mates; no per-base / windowed output",
	 "Quick per-chrom summary incl. base/map quality; one-liner sanity check"],

	["samtools idxstats (htslib / pysam)",
	 "**Python / pysam** in the uv venv — no Docker<br>`pysam.idxstats(CRAM)`<br>_(or samtools CLI)_",
	 "Not needed",
	 "Instant on BAM; ≈35 s on CRAM (.crai lacks counts → scans)",
	 "refname, seqlen, mapped, unmapped (read counts)",
	 "Per-chromosome only",
	 "Counts, not depth (÷length ×read_len for est.); NO flag filtering (incl. dups/secondary → high); not instant on CRAM; multimapping inflates chrY/repeats",
	 "Ultra-fast per-chrom read counts (BAM); rough CNV via reads/base"],

	["bamtocov / bamtocounts (Nim / hts-nim)",
	 "**Docker biocontainer** `quay.io/biocontainers/bamtocov:2.8.0--h1104d80_0`<br>`bamtocov [--physical] BAM` (BAM only); `bamtocounts --fasta REF TARGET CRAM`<br>_(also: bioconda, Nim source)_",
	 "bamtocov ✗ (BAM only); bamtocounts `--fasta` ✓",
	 "Streams whole file (no indexed seek)",
	 "bamtocov: per-base / WIG / quantized, **physical coverage**, stranded; bamtocounts: per-target counts",
	 "bamtocov `--wig SPAN` / `--regions`; bamtocounts needs target BED/GFF",
	 "Main binary can't read CRAM (convert first); bamtocounts reads entire file; `-r` means different things across binaries",
	 "Physical coverage (SV / structural QC), stranded / quantized depth classes, per-feature counts"],

	["megadepth (C++; recount3)",
	 "**Docker biocontainer** `quay.io/biocontainers/megadepth:1.2.0--hff880f7_4`<br>`megadepth CRAM --fasta REF [--annotation BED --op mean] [--bigwig]`<br>_(also: bioconda, precompiled binary, pip `megadepth`)_",
	 "`--fasta`",
	 "Very fast",
	 "BigWig; per-annotation / region summaries; AUC; per-base coverage",
	 "Via `--annotation` BED (no fixed-window flag)",
	 "No built-in fixed-window binning; oriented to BigWig / feature quantification",
	 "BigWig tracks for browsers; gene/exon quantification (RNA-seq style); AUC normalization"],

	["goleft indexcov (Go)",
	 "**Docker biocontainer** `quay.io/biocontainers/goleft:0.2.6--he881be0_1`<br>`goleft indexcov -d OUT --fai REF.fai --extranormalize CRAM.crai`<br>_(also: bioconda, single-file Go binary)_",
	 "Needs `--fai`; reads `.crai`",
	 "Fastest — ≈1.8 s, index-only",
	 "Interactive HTML (overview + per-chrom); `.bed.gz` ~16 kb-bin scaled cov (~1.0 = diploid); `.ped` inferred sex + X/Y CN; `.roc`",
	 "Fixed ~16,384 bp bins from the index",
	 "Coarse (~16 kb); estimate, not exact depth; relative (scaled) not absolute; multimapping inflates repeats; needs `.fai` for CRAM",
	 "Near-instant genome-wide CNV / aneuploidy scan; screening many samples; QC first-pass"],
]

# Specialized / analysis-specific coverage approaches examined (not interchangeable, not run)
SPECIAL_COLUMNS = ["Tool", "What 'coverage' means here", "Method", "Use case"]
SPECIAL_ROWS = [
	["bamsnap",
	 "Per-base depth + A/C/G/T composition over a tiny region",
	 "Per-read iteration into a Python dict over a small window (for rendering)",
	 "Variant snapshot images (visualization), not genome-wide"],
	["EchoSV `get_coverage_bam`",
	 "Uniqueness-aware per-base depth, SV-oriented",
	 "numpy per-base array; MAPQ≥1; CIGAR-aware (small del ≤20 bp counted, big del/N = gap); RLE region output; also builds read→ref segment chains",
	 "Structural-variant coverage from haplotype/contig-to-reference alignments; exact breakpoints"],
	["PaSDqc",
	 "Frequency content of the depth signal (magnitude discarded)",
	 "Power spectral density (Lomb–Scargle + Welch) of mean-normalized depth at uniquely-mappable positions",
	 "scWGS amplification-uniformity QC; amplicon-size distribution; read-depth variance"],
]


def _md_table(columns, rows):
	out = ["| " + " | ".join(columns) + " |",
	       "| " + " | ".join(["---"] * len(columns)) + " |"]
	for r in rows:
		# escape pipes inside cells
		cells = [c.replace("|", "\\|") for c in r]
		out.append("| " + " | ".join(cells) + " |")
	return "\n".join(out)


def write_markdown(path):
	md = []
	md.append("# Coverage tool comparison — ParkLab ChromCov\n")
	md.append("Comparison of read-coverage software evaluated against the COLO829T CRAM "
	          "(GRCh38, ~15× autosomal). All practical tools were run via bioconda "
	          "biocontainers (`quay.io/biocontainers/...`).\n")
	md.append("## Practical coverage tools\n")
	md.append(_md_table(COLUMNS, ROWS))
	md.append("\n## Specialized / analysis-specific approaches\n")
	md.append("Examined for how they compute coverage, but not interchangeable depth tools:\n")
	md.append(_md_table(SPECIAL_COLUMNS, SPECIAL_ROWS))
	md.append("\n## Cross-cutting notes\n")
	md.append("- **Delivery:** pysam tools (samtools coverage / idxstats) install into the uv "
	          "venv (the pysam wheel bundles htslib) — no Docker. The compiled tools (mosdepth, "
	          "PanDepth, bamtocov, megadepth, goleft) run via Docker biocontainers here, since "
	          "they aren't pip-installable and GitHub binaries are Linux-only on this arm64 Mac.\n"
	          "- **Overlap convention:** mosdepth is mate-overlap aware; samtools/PanDepth "
	          "count both mates → ~0.93× constant offset. Cancels in CNV ratios; matters for "
	          "absolute depth / variant support.\n"
	          "- **Index-only speed:** truly instant per-chrom counts need a BAM `.bai`; a CRAM "
	          "`.crai` lacks per-reference counts, so `idxstats` scans (~35 s). `goleft indexcov` "
	          "is the fast index-based route (binned coverage from `.crai`).\n"
	          "- **Filtering:** depth tools exclude flag 1796 (dup/secondary/supplementary/QC-fail); "
	          "`idxstats` counts everything → reads high, especially on repetitive chroms (chrY).\n"
	          "- **For CNV/relative work** any tool agrees; **for absolute depth** mind overlap + "
	          "filtering conventions.\n")
	path.write_text("\n".join(md))
	return path


def _inline_md_to_html(text):
	text = html.escape(text)
	text = text.replace("&lt;br&gt;", "<br>")          # keep explicit line breaks
	text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
	text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
	text = re.sub(r"_([^_]+)_", r"<em>\1</em>", text)
	return text


def _html_table(columns, rows):
	head = "".join(f"<th>{_inline_md_to_html(c)}</th>" for c in columns)
	body = ""
	for r in rows:
		tds = "".join(f"<td>{_inline_md_to_html(c)}</td>" for c in r)
		body += f"<tr>{tds}</tr>\n"
	return f"<table>\n<thead><tr>{head}</tr></thead>\n<tbody>\n{body}</tbody>\n</table>"


def write_html(path):
	css = """
	body { font-family: -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
	       margin: 2rem; color: #1a1a1a; }
	h1 { margin-bottom: .2rem; } h2 { margin-top: 2rem; }
	p.lead { color: #555; }
	table { border-collapse: collapse; width: 100%; font-size: 13px; margin-top: .5rem; }
	th, td { border: 1px solid #d6dbe0; padding: 8px 10px; vertical-align: top; text-align: left; }
	thead th { position: sticky; top: 0; background: #2c3e50; color: #fff; }
	tbody tr:nth-child(even) { background: #f5f7fa; }
	td:first-child, th:first-child { font-weight: 600; white-space: nowrap; }
	code { background: #eef1f4; padding: 1px 4px; border-radius: 3px;
	       font-family: SFMono-Regular, Menlo, Consolas, monospace; font-size: 12px; }
	ul { font-size: 14px; } li { margin: .25rem 0; }
	"""
	parts = [f"<!doctype html><html><head><meta charset='utf-8'>",
	         "<title>Coverage tool comparison — ParkLab ChromCov</title>",
	         f"<style>{css}</style></head><body>",
	         "<h1>Coverage tool comparison — ParkLab ChromCov</h1>",
	         "<p class='lead'>Read-coverage software evaluated against the COLO829T CRAM "
	         "(GRCh38, ~15× autosomal). Practical tools run via bioconda biocontainers.</p>",
	         "<h2>Practical coverage tools</h2>", _html_table(COLUMNS, ROWS),
	         "<h2>Specialized / analysis-specific approaches</h2>",
	         "<p class='lead'>How they compute coverage — not interchangeable depth tools.</p>",
	         _html_table(SPECIAL_COLUMNS, SPECIAL_ROWS),
	         "<h2>Cross-cutting notes</h2><ul>",
	         "<li><strong>Delivery:</strong> pysam tools (samtools coverage / idxstats) install "
	         "into the uv venv (the pysam wheel bundles htslib) — no Docker. The compiled tools "
	         "(mosdepth, PanDepth, bamtocov, megadepth, goleft) run via Docker biocontainers here, "
	         "since they aren't pip-installable and GitHub binaries are Linux-only on this arm64 "
	         "Mac.</li>",
	         "<li><strong>Overlap convention:</strong> mosdepth is mate-overlap aware; "
	         "samtools/PanDepth count both mates → ~0.93× constant offset. Cancels in CNV ratios; "
	         "matters for absolute depth / variant support.</li>",
	         "<li><strong>Index-only speed:</strong> instant per-chrom counts need a BAM "
	         "<code>.bai</code>; a CRAM <code>.crai</code> lacks per-reference counts so "
	         "<code>idxstats</code> scans (~35 s). <code>goleft indexcov</code> is the fast "
	         "index-based route.</li>",
	         "<li><strong>Filtering:</strong> depth tools exclude flag 1796 "
	         "(dup/secondary/supplementary/QC-fail); <code>idxstats</code> counts everything → "
	         "reads high, especially on repetitive chroms (chrY).</li>",
	         "<li><strong>Rule of thumb:</strong> for CNV/relative work any tool agrees; for "
	         "absolute depth, mind overlap + filtering conventions.</li>",
	         "</ul></body></html>"]
	path.write_text("".join(parts))
	return path


if __name__ == "__main__":
	md = write_markdown(DOCS / "coverage_tools_comparison.md")
	ht = write_html(DOCS / "coverage_tools_comparison.html")
	print(f"wrote {md}")
	print(f"wrote {ht}")
