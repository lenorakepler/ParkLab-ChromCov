"""
chromcov — per-chromosome average coverage from a CRAM.

Public API:
  Config                           configuration (Pydantic, validated at the edge)
  pipeline.run(cfg, depth=...)     the single orchestrator -> RunResult
  Depth / Source                   run parameters (MEAN vs FULL; CRAM vs tracks)
  RunResult                        the accumulator a run folds into
  ChromCoverage                    one normalized per-chromosome row
  ChromDepth / ChromStats / DepthHistogram   the reduction classes
  Strata                           callability categories (masking mechanism)
  QCThresholds                     tunable abnormality-flag thresholds

per-base track I/O lives in `chromcov.io.track`; mosdepth is an optional
cross-check add-on (scripts/mosdepth_coverage.py).
"""
from .config.schema import Config, QCThresholds
from .categories import Strata
from .pipeline import Depth, Source, run
from .present.frames import ChromCoverage
from .reduce import ChromDepth, ChromStats, DepthHistogram
from .result import RunResult

__all__ = [
    "Config",
    "QCThresholds",
    "Depth",
    "Source",
    "run",
    "RunResult",
    "ChromCoverage",
    "ChromDepth",
    "ChromStats",
    "DepthHistogram",
    "Strata",
]
