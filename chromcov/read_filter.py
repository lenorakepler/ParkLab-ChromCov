"""
Configurable read filtering for coverage calculation.

This is a standalone variant of `get_cov.py` where the hardcoded read filters in
`does_fail` are replaced by a user-configurable `ReadFilter`. Filters can be set from a CLI
or config file and are translated *once* into an efficient representation, so the per-read
cost in the hot fetch loop stays O(1) no matter how many filters are enabled.

Key idea: `read.is_unmapped`, `read.is_secondary`, etc. are just convenience accessors over
bits of the single integer `read.flag`. So any number of flag filters collapse into one
integer bitmask -- exactly how samtools/mosdepth do it. Flag semantics mirror samtools:
  -f  include_flags     : keep only reads with ALL of these bits set
  -F  exclude_flags     : drop reads with ANY of these bits set
  -G  exclude_all_flags : drop reads only if ALL of these bits are set
plus a MAPQ threshold (-Q).
"""
from pathlib import Path
from dataclasses import dataclass
import pysam
import numpy as np

# SAM flag name -> bit, so config/CLI can use friendly names while internals use integers.
SAM_FLAGS = {
    "paired":        0x1,
    "proper_pair":   0x2,
    "unmapped":      0x4,
    "mate_unmapped": 0x8,
    "reverse":       0x10,
    "mate_reverse":  0x20,
    "read1":         0x40,
    "read2":         0x80,
    "secondary":     0x100,
    "qcfail":        0x200,
    "duplicate":     0x400,
    "supplementary": 0x800,
}

# unmapped | secondary | duplicate | supplementary -- reproduces get_cov.does_fail defaults.
DEFAULT_EXCLUDE = SAM_FLAGS["unmapped"] | SAM_FLAGS["secondary"] | SAM_FLAGS["duplicate"] | SAM_FLAGS["supplementary"]


def _to_mask(flags):
    """Normalize a flag spec to an int bitmask.

    Accepts an int mask (returned as-is), an iterable of flag names/ints (OR-ed together),
    or None/empty (-> 0). Lets config pass either a raw integer (samtools parity) or names
    like ["unmapped", "secondary"].
    """
    if flags is None:
        return 0
    if isinstance(flags, int):
        return flags
    mask = 0
    for f in flags:
        mask |= f if isinstance(f, int) else SAM_FLAGS[f]
    return mask


@dataclass
class ReadFilter:
    """Read-level filter built once from config, then applied cheaply per read.

    Flag args accept either an int mask or an iterable of names/ints (see `_to_mask`).
    """
    include_flags: int = 0                    # -f: require ALL these bits
    exclude_flags: int = DEFAULT_EXCLUDE      # -F: exclude if ANY set (default = old behavior)
    exclude_all_flags: int = 0                # -G: exclude only if ALL set
    min_mapping_quality: int = 0              # -Q

    def __post_init__(self):
        self.include_flags     = _to_mask(self.include_flags)
        self.exclude_flags     = _to_mask(self.exclude_flags)
        self.exclude_all_flags = _to_mask(self.exclude_all_flags)

    def fails(self, read):
        flag = read.flag
        if flag & self.exclude_flags:                                  # -F: any excluded bit set
            return True
        if (flag & self.include_flags) != self.include_flags:          # -f: missing a required bit
            return True
        if self.exclude_all_flags and (flag & self.exclude_all_flags) == self.exclude_all_flags:  # -G
            return True
        if read.mapping_quality < self.min_mapping_quality:            # -Q
            return True
        return False


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
        base_depth = np.cumsum(tally[:-1])   # per-base depth, zeros included for free
        total_depth = np.sum(base_depth)

    else:
        base_depth = None
        total_depth = tally

    chrom_cov = total_depth / chrom_length

    return base_depth, total_depth, chrom_cov


if __name__ == "__main__":
    import timeit

    data_dir = Path().resolve() / "data"

    cram_file = str(data_dir / "COLO829T_TEST.cram")
    cram_ref = str(data_dir / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa")
    cram_index = str(data_dir / "COLO829T_TEST.cram.crai")

    cram = pysam.AlignmentFile(cram_file, "rc", reference_filename=cram_ref, index_filename=cram_index)

    chroms = ["chrUn_KI270382v1", "chr21"]

    # Default filter reproduces the old hardcoded does_fail behavior.
    read_filter = ReadFilter()

    def per_base():
        for chrom in chroms:
            base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, read_filter, per_base=True)
            print(f"Per base, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")

    def aggregated():
        for chrom in chroms:
            base_depth, total_depth, chrom_cov = calc_cov(cram, chrom, read_filter, per_base=False)
            print(f"Aggregated, {chrom}: {base_depth}, {total_depth}, {chrom_cov}")

    base_time = timeit.timeit(per_base, number=1)
    print(f"--> Time for per base: {base_time:.3f} seconds\n")

    agg_time = timeit.timeit(aggregated, number=1)

    print(f"--> Time for aggregated: {agg_time:.3f} seconds")
