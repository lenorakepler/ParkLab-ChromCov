"""
Coverage entry point: preflight, then the hand-rolled pysam calculator.

There is a single in-tool calculator (the native `calc_cov`). mosdepth is kept as
an OPTIONAL cross-check add-on that emits the same output format
(`scripts/mosdepth_coverage.py`), not a second backend baked into the dispatch
path -- so the core stays free of subprocess/registry machinery.
"""
from __future__ import annotations

import pysam

from . import validate
from .config import CoverageConfig
from .read_filter import ReadFilter, calc_cov
from .result import ChromCoverage


def run_coverage(config: CoverageConfig, skip_preflight: bool = False) -> list[ChromCoverage]:
    """Per-chromosome coverage rows via the native calculator. Kept free of
    argparse/I-O so tests can call it directly."""
    if not skip_preflight:
        validate.preflight(config)   # sorted / indexed / reference-M5; raises on failure

    rf = ReadFilter(
        include_flags=config.include_flags,
        exclude_flags=config.exclude_flags,
        exclude_all_flags=config.exclude_all_flags,
        min_mapping_quality=config.min_mapping_quality,
    )
    cram = pysam.AlignmentFile(
        str(config.cram), "rc",
        reference_filename=str(config.reference),
        index_filename=str(config.index),
    )
    try:
        chroms = config.chroms if config.chroms is not None else cram.references
        rows: list[ChromCoverage] = []
        for chrom in chroms:
            _, total_depth, _ = calc_cov(cram, chrom, rf, per_base=config.per_base)
            rows.append(ChromCoverage(
                chrom=chrom,
                length=cram.get_reference_length(chrom),
                bases=int(total_depth),
            ))
    finally:
        cram.close()
    return rows
