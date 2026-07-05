"""
chromcov — per-chromosome average coverage from a CRAM.

Public API:
  CoverageConfig / AnalysisConfig / RunConfig   configuration (Pydantic, validated)
  run_coverage(config)             per-chromosome coverage -> list[ChromCoverage]
  ChromCoverage                    one normalized per-chromosome row
  CoverageAnalysis                 the full QC pipeline (stats/windows/strata/plots)
  ChromDepth / Strata              the reduction + callability classes
  PerBaseStore / RunStore          per-base track store + coverage-table archive

(mosdepth is an optional cross-check add-on: scripts/mosdepth_coverage.py.)
"""
from .analysis import ChromDepth, ChromStats, DepthHistogram
from .config import AnalysisConfig, CoverageConfig, RunConfig
from .dispatch import run_coverage
from .output import RunStore
from .perbase import PerBaseStore
from .pipeline import CoverageAnalysis
from .qc import QCThresholds
from .result import ChromCoverage
from .strata import Strata

__all__ = [
    "CoverageConfig",
    "AnalysisConfig",
    "RunConfig",
    "run_coverage",
    "ChromCoverage",
    "CoverageAnalysis",
    "ChromDepth",
    "ChromStats",
    "DepthHistogram",
    "Strata",
    "RunStore",
    "PerBaseStore",
    "QCThresholds",
]
