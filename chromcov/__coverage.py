from pathlib import Path
import pandas as pd
import pysam

DATA_DIR = Path().resolve() / "data"
OUT_DIR = Path().resolve() / "out"
OUT_DIR.mkdir(exist_ok=True)

cram_file = str(DATA_DIR / "COLO829T_TEST.cram")
cram_ref = str(DATA_DIR / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")

main_chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Assumed Illumina read length, used to turn idxstats read counts into a depth
# estimate (depth ~= reads * READ_LEN / chrom_length).
READ_LEN = 151

# Columns emitted by `samtools coverage` that are numeric
NUM_COLS = ["startpos", "endpos", "numreads", "covbases",
            "coverage", "meandepth", "meanbaseq", "meanmapq"]


def coverage_per_chrom(cram=cram_file, ref=cram_ref, chroms=main_chroms):
	"""Per-chromosome coverage via `samtools coverage`, dispatched through pysam.

	Loops one indexed region per chromosome so only the main chromosomes are
	read (skipping the ~170 unplaced/alt contigs). Returns a tidy DataFrame.
	"""
	rows = []
	for chrom in chroms:
		# pysam.coverage returns the samtools coverage table as a string:
		# a "#"-prefixed header line plus one data row for the region.
		out = pysam.coverage("--reference", ref, "-r", chrom, cram)
		header_line, data_line = out.strip().split("\n")
		header = header_line.lstrip("#").split("\t")
		values = data_line.split("\t")
		rows.append(dict(zip(header, values)))
		print(f"{chrom}: done", flush=True)

	df = pd.DataFrame(rows)
	for c in NUM_COLS:
		df[c] = pd.to_numeric(df[c])
	return df


def build_comparison():
	"""Merge per-chromosome mean depth from every method that has produced
	output, into one comparison table. Each source is optional: a method is
	skipped if its file isn't present yet.
	"""
	cov = pd.DataFrame({"chrom": main_chroms}).set_index("chrom")

	# samtools coverage (this module)
	f = OUT_DIR / "coverage_samtools.csv"
	if f.exists():
		s = pd.read_csv(f).set_index("rname")["meandepth"]
		cov["samtools"] = s.reindex(main_chroms)

	# PanDepth (.chr.stat.gz: #Chr, Length, CoveredSite, TotalDepth, Coverage(%), MeanDepth)
	f = DATA_DIR / "pandepth_chrom.chr.stat.gz"
	if f.exists():
		p = pd.read_csv(f, sep="\t").rename(columns={"#Chr": "chrom"}).set_index("chrom")["MeanDepth"]
		cov["pandepth"] = p.reindex(main_chroms)

	# mosdepth (summary.txt: chrom, length, bases, mean, min, max)
	f = DATA_DIR / "testmosslow.mosdepth.summary.txt"
	if f.exists():
		m = pd.read_csv(f, sep="\t").set_index("chrom")["mean"]
		cov["mosdepth"] = m.reindex(main_chroms)

	# samtools idxstats (refname, seqlen, mapped, unmapped). Convert read counts
	# to a depth estimate: depth ~= mapped * READ_LEN / seqlen, and a CNV ratio
	# from baseline-normalized reads-per-base. NB: idxstats counts ALL mapped
	# alignments (incl. duplicates/secondary), so idx_depth reads slightly high;
	# the ratio is unaffected. On CRAM this is not instant (the .crai lacks counts).
	f = OUT_DIR / "idxstats.tsv"
	if f.exists():
		ix = pd.read_csv(f, sep="\t", names=["chrom", "seqlen", "mapped", "unmapped"])
		ix = ix.set_index("chrom").reindex(main_chroms)
		cov["idx_depth"] = (ix["mapped"] * READ_LEN / ix["seqlen"]).round(2)
		autosomes = [f"chr{i}" for i in range(1, 23)]
		base_rpb = (ix.loc[autosomes, "mapped"] / ix.loc[autosomes, "seqlen"]).mean()
		cov["idx_ratio"] = ((ix["mapped"] / ix["seqlen"]) / base_rpb).round(3)

	# goleft indexcov (bed.gz: #chrom, start, end, <sample>) -> scaled coverage,
	# already normalized so ~1.0 = diploid. Mean over the ~16 kb bins per chrom.
	f = OUT_DIR / "indexcov" / "indexcov-indexcov.bed.gz"
	if f.exists():
		ic = pd.read_csv(f, sep="\t")
		ic.columns = ["chrom", "start", "end", "scaled"]
		cov["indexcov"] = ic.groupby("chrom")["scaled"].mean().reindex(main_chroms).round(3)

	# CNV-style ratio from the autosomal length-weighted baseline (uses pandepth if present, else first available)
	depth_col = next((c for c in ["pandepth", "mosdepth", "samtools"] if c in cov), None)
	if depth_col:
		autosomes = [f"chr{i}" for i in range(1, 23)]
		baseline = cov.loc[autosomes, depth_col].mean()
		cov["ratio"] = cov[depth_col] / baseline
		cov["approx_cn"] = 2 * cov["ratio"]

	out_path = OUT_DIR / "coverage_comparison.csv"
	cov.to_csv(out_path)
	print(cov.round(2).to_string())
	print(f"wrote {out_path}")
	return cov


if __name__ == "__main__":
	import sys
	if "--merge-only" not in sys.argv:
		df = coverage_per_chrom()
		out_path = OUT_DIR / "coverage_samtools.csv"
		df.to_csv(out_path, index=False)
		print(df[["rname", "numreads", "covbases", "coverage", "meandepth"]].to_string(index=False))
		print(f"wrote {out_path}")
	build_comparison()
