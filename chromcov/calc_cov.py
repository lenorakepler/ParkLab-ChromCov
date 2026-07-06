"""
Set up filtering out reads:
	class: ReadFilter
	_to_mask()

Main coverage calculation:
	calc_cov()
"""
from pathlib import Path

import numpy as np
import pysam

def calc_cov(cram, chrom, read_filter, per_base=False):
    chrom_length = cram.get_reference_length(chrom)

    if per_base:
        # running finite difference array,
        # +1 to account for last pos off chromosome
        tally = np.zeros(chrom_length + 1, dtype=np.int32)
    else:
        # running total of bases
        tally = 0

    # Iterate over each read in the alignment
    for read in cram.fetch(chrom):

        # ignore reads that do not pass filter criteria
        if read_filter.fails(read):
            continue

        # get_blocks iterates through aligned gapless blocks
        #   0-based half-open [start, end)
        #   defined by CIGAR ops of M, =, or X (aln match, seq match, seq mismatch)
        #   no I, D, N, S, or H (ins, del, ref skip, soft clip, hard clip)
        #   i.e. they match both sequence and reference
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


if __name__ == "__main__":
    import timeit

    from chromcov.config import ReadFilter

    data_dir = Path().resolve() / "data"
    cram_file = str(data_dir / "COLO829T_TEST.cram")
    cram_ref = str(data_dir / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
    cram_index = str(data_dir / "COLO829T_TEST.cram.crai")

    cram = pysam.AlignmentFile(cram_file, "rc", reference_filename=cram_ref, index_filename=cram_index)
    chroms = ["chrUn_KI270382v1", "chr21"]
    read_filter = ReadFilter()

    def per_base():
        for chrom in chroms:
            base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, read_filter, per_base=True)
            print(f"Per base, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")

    def aggregated():
        for chrom in chroms:
            base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, read_filter, per_base=False)
            print(f"Aggregated, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")

    print(f"--> Time for per base: {timeit.timeit(per_base, number=1):.3f} seconds\n")
    print(f"--> Time for aggregated: {timeit.timeit(aggregated, number=1):.3f} seconds")
