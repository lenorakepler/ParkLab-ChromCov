"""
chromcov — per-chromosome average coverage from a CRAM.

Public API:
  CoverageConfig / AnalysisConfig  configuration (Pydantic, validated at the boundary)
  run_coverage(config)             dispatch to a backend -> list[ChromCoverage]
  CoverageBackend / get_backend    the backend ABC + registry
  ChromCoverage                    one normalized per-chromosome row
  CoverageAnalysis                 the full QC pipeline (stats/windows/strata/plots)
  ChromDepth / Strata / RunStore   the reduction, callability, and archival classes
"""
from .analysis import ChromDepth, ChromStats
from .backends import CoverageBackend, MosdepthBackend, NativeBackend, get_backend
from .config import AnalysisConfig, CoverageConfig
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
    "run_coverage",
    "CoverageBackend",
    "NativeBackend",
    "MosdepthBackend",
    "get_backend",
    "ChromCoverage",
    "CoverageAnalysis",
    "ChromDepth",
    "ChromStats",
    "Strata",
    "RunStore",
    "PerBaseStore",
    "QCThresholds",
]
