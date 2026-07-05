"""
Backend dispatch. Deliberately trivial: pick a callable, both share the
`run(config) -> list[ChromCoverage]` contract. Keep this free of argparse/I-O so
tests can call `run_coverage(config)` directly and so "run both, diff the means"
is a 3-line cross-validation test.
"""
from config import CoverageConfig
from result import ChromCoverage

import native
import mosdepth

BACKENDS = {
    "native": native.run,
    "mosdepth": mosdepth.run,
}


def run_coverage(config: CoverageConfig) -> list[ChromCoverage]:
    try:
        backend = BACKENDS[config.backend]
    except KeyError:
        raise ValueError(f"unknown backend {config.backend!r}; choose from {sorted(BACKENDS)}")
    return backend(config)
