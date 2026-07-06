"""
chromcov — per-chromosome average coverage from a CRAM.

Public API:
  Config                           configuration (Pydantic, validated at the edge)
  run_coverage(config)             per-chromosome MEAN coverage -> list[ChromCoverage]
  ChromCoverage                    one normalized per-chromosome row
  QCReport                         the --full driver (compute bedgraphs -> reduce)
  ChromDepth / Strata              the reduction + callability classes

per-base bedgraph I/O lives in `chromcov.perbase`; mosdepth is an optional
cross-check add-on (scripts/mosdepth_coverage.py).
"""
from .config import Config
from .coverage import run_coverage
from .depth import ChromDepth, ChromStats, DepthHistogram
from .config import QCThresholds
from .qc_report import QCReport
from .report import ChromCoverage
from .strata import Strata

__all__ = [
    "Config",
    "run_coverage",
    "ChromCoverage",
    "QCReport",
    "ChromDepth",
    "ChromStats",
    "DepthHistogram",
    "Strata",
    "QCThresholds",
]
