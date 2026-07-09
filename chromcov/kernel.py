"""
The coverage kernel: `calc_cov` tallies aligned bases over one contig.

Read filtering (`ReadFilter`, `to_mask`) lives in `chromcov.filtering`; this
module is only the depth computation. `per_base=False` returns just the base
total (the mean path); `per_base=True` also returns the exact per-base depth
vector (the --full path) -- the two totals agree because the mean is the same
number computed without materializing the vector.
"""
import numpy as np


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
