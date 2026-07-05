"""
Native backend: thin adapter over the hand-rolled pysam calculator.

Translates the shared CoverageConfig into a `chromcov.read_filter.ReadFilter`
once, opens the CRAM once, and loops chromosomes -- emitting the same
`list[ChromCoverage]` the mosdepth backend does. All the actual coverage math
stays in the package (`calc_cov`); this file only bridges config <-> that API.
"""
from __future__ import annotations
import sys
from pathlib import Path

# Sketch lives outside the package and uses flat imports for its siblings, so put
# the repo root on the path to reach `chromcov`. (In the real package this file
# would live under chromcov/backends/ and just do a normal relative import.)
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import pysam
from chromcov.read_filter import ReadFilter, calc_cov

from config import CoverageConfig
from result import ChromCoverage


def run(config: CoverageConfig) -> list[ChromCoverage]:
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

    chroms = config.chroms if config.chroms is not None else cram.references

    rows: list[ChromCoverage] = []
    for chrom in chroms:
        base_depth, total_depth, _ = calc_cov(cram, chrom, rf, per_base=config.per_base)
        rows.append(
            ChromCoverage(
                chrom=chrom,
                length=cram.get_reference_length(chrom),
                bases=int(total_depth),
                backend="native",
                base_depth=base_depth,
            )
        )
    return rows
