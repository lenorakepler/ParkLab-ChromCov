"""
Backend dispatch. Deliberately trivial: resolve a backend from the registry,
run preflight (unless skipped), delegate. Kept free of argparse/I-O so tests can
call `run_coverage(config)` directly and so "run both, diff the means" is a
3-line cross-validation.
"""
from __future__ import annotations

from . import validate
from .backends import get_backend
from .config import CoverageConfig
from .result import ChromCoverage


def run_coverage(config: CoverageConfig, skip_preflight: bool = False) -> list[ChromCoverage]:
    backend = get_backend(config.backend)
    if not skip_preflight:
        validate.preflight(config)   # sorted / indexed / reference-M5; raises on failure
    return backend.run(config)
