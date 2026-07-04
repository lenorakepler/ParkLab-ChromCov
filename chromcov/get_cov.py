"""
mosdepth-like coverage calculation using event accumulation / finite difference array
(+/- 1 at start/stop of alignment blocks --> cumulative sum = per-base coverage without incrementing each position)
"""
from pathlib import Path
import pysam
import numpy as np

def does_fail(read, min_map_qual, **fail_criteria):
	"""
	read.is_unmapped = True
	read.is_secondary = True
	read.is_supplementary = True
	read.is_duplicate = True
	read.mapping_quality < min_map_qual
	read.qual < min_base_qual

	read.flag in [filter_flags]
	"""
	if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_duplicate:
		return True

	if read.mapping_quality < min_map_qual:
		return True

	return False

def calc_cov(cram, chrom, min_map_qual=0, per_base=False):
	"""
	mosdepth-like coverage calculation using event accumulation / finite difference array
	(+/- 1 at start/stop of alignment blocks --> cumulative sum = per-base coverage without incrementing each position)
	"""
	chrom_length = cram.get_reference_length(chrom)
	
	if per_base:
		# running finite difference array, 
		# +1 to account for last pos off chromosome
		tally = np.zeros(chrom_length + 1, dtype=np.int32)

	else:
		# running total of bases
		tally = 0

	for read in cram.fetch(chrom):
		if does_fail(read, min_map_qual):
			continue

		# get_blocks iterates through aligned gapless blocks 
		#   0-based half-open [start, end)
		#   defined by CIGAR ops of M, =, or X (aln match, seq match, seq mismatch)
		#   no I, D, N, S, or H (ins, del, ref skip, soft clip, hard clip)
		#   i.e. they match both sequence and reference
		#
		# NOTE: unlike pileup(), this DOES NOT ACCOUNT FOR OVERLAPING MATE PAIRS
		for start, end in read.get_blocks():
			if per_base:
				tally[start] += 1
				tally[end]   -= 1
			else:
				tally += end - start

	if per_base:
		# Per-base depth is cumulative sum of finite diff array
		# (leaving off added index for first off-chrom base)
		base_depth = np.cumsum(tally[:-1])
		total_depth = np.sum(base_depth)

	else:
		base_depth = None
		total_depth = tally
	
	chrom_cov = total_depth / chrom_length
	return base_depth, total_depth, chrom_cov

def breadth_at_depth(base_depth, bins=10):
	# Take the cumulative sum in reverse order
	# (start at highest depth)
	hist, bin_edges = np.histogram(base_depth, bins=bins)
	cum_breadth_rev = np.cumsum(hist[::-1])

	# Reverse
	cum_breadth = cum_breadth_rev[::-1]
	cum_breadth_pct = cum_breadth / cum_breadth[0]

	return breadth_hist, breadth_pcts, bin_edges
	
if __name__ == "__main__":
	import timeit

	data_dir = Path().resolve() / "data"

	cram_file = str(data_dir / "COLO829T_TEST.cram")
	cram_ref = str(data_dir / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
	cram_index = str(data_dir / "COLO829T_TEST.cram.crai")

	cram = pysam.AlignmentFile(cram_file, "rc", reference_filename=cram_ref, index_filename=cram_index)

	chroms = ["chrUn_KI270382v1", "chr21"]

	def per_base():
		for chrom in chroms:
			base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, min_map_qual=0, per_base=True)
			print(f"Per base, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")
			breakpoint()

	def aggregated():
		for chrom in chroms:
			base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, min_map_qual=0, per_base=False)
			print(f"Aggregated, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")

	base_time = timeit.timeit(per_base, number=1)
	print(f"--> Time for per base: {base_time:.3f} seconds\n")

	agg_time = timeit.timeit(aggregated, number=1)
	print(f"--> Time for aggregated: {agg_time:.3f} seconds")

