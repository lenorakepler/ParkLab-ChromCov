"""
Validate input files and run coverage calculation
"""
from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
import pysam

from . import validate
from .calc_cov import calc_cov
from .config import Config, ReadFilter
from .report import ChromCoverage

def _read_filter(cfg: Config) -> ReadFilter:
    return ReadFilter(
        include_flags=cfg.include_flags,
        exclude_flags=cfg.exclude_flags,
        exclude_all_flags=cfg.exclude_all_flags,
        min_mapping_quality=cfg.min_mapping_quality,
    )

def _mean_one(cfg: Config, chrom: str) -> ChromCoverage:
    cram = pysam.AlignmentFile(
        str(cfg.cram), "rc",
        reference_filename=str(cfg.reference), index_filename=str(cfg.index))
    try:
        _, total_depth, _ = calc_cov(cram, chrom, _read_filter(cfg), per_base=False)
        length = cram.get_reference_length(chrom)
    finally:
        cram.close()
    return ChromCoverage(chrom=chrom, length=length, bases=int(total_depth))


def _mean_worker(args) -> ChromCoverage:
    cfg, chrom = args
    return _mean_one(cfg, chrom)

def run_coverage(config: Config, skip_preflight: bool = False, jobs: int = 1) -> list[ChromCoverage]:
    """
    Per-chromosome mean coverage. 
    Kept free of argparse/IO so tests/workflows can call it directly.
    `jobs > 1` computes contigs in parallel.
    """
    if not skip_preflight:
        validate.preflight(config)   # sorted / indexed / reference-M5; raises on failure

    cram = pysam.AlignmentFile(
        str(config.cram), "rc",
        reference_filename=str(config.reference), index_filename=str(config.index))
    try:
        chroms = config.select_contigs(cram.references)
    finally:
        cram.close()

    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            return list(ex.map(_mean_worker, [(config, c) for c in chroms]))

    # Serial: reuse a single open handle across contigs (cheaper for many contigs).
    rf = _read_filter(config)
    cram = pysam.AlignmentFile(
        str(config.cram), "rc",
        reference_filename=str(config.reference), index_filename=str(config.index))
    try:
        rows: list[ChromCoverage] = []
        for chrom in chroms:
            _, total_depth, _ = calc_cov(cram, chrom, rf, per_base=False)
            rows.append(ChromCoverage(
                chrom=chrom, length=cram.get_reference_length(chrom), bases=int(total_depth)))
    finally:
        cram.close()
    return rows
