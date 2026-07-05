"""
Cross-validation: the hand-rolled native calculator vs. mosdepth.

The two backends should agree because they count the same thing -- aligned bases
over M/=/X CIGAR ops (D/N gaps excluded), no overlapping-mate correction, and
`mean = bases / length`. With `exclude_flags` pinned to the same mask in
CoverageConfig, any disagreement is a real bug in one of them, which is exactly
what makes running both worthwhile.

Run (pytest isn't a project dependency):

    uv run --with pytest pytest dev/dispatch-sketch/test_backends.py -v

Skips cleanly when mosdepth or the data/ CRAM aren't present, so it's safe in CI
on a machine without either.
"""
from __future__ import annotations
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))

from config import CoverageConfig
import dispatch

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATA = _REPO_ROOT / "data"
_CRAM = _DATA / "COLO829T_TEST.cram"
_REF = _DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# A small contig keeps the mosdepth subprocess fast; chr21 is the smallest
# standard autosome and is CNV-interesting in COLO829T.
_CONTIG = "chr21"

requires_mosdepth = pytest.mark.skipif(
    shutil.which("mosdepth") is None, reason="mosdepth not on PATH"
)
requires_data = pytest.mark.skipif(
    not (_CRAM.exists() and _REF.exists()), reason="data/ CRAM or reference missing"
)


def _config(backend: str) -> CoverageConfig:
    return CoverageConfig(
        cram=_CRAM,
        reference=_REF,
        backend=backend,
        chroms=(_CONTIG,),
    )


@requires_data
@requires_mosdepth
def test_native_matches_mosdepth():
    (native,) = dispatch.run_coverage(_config("native"))
    (mos,) = dispatch.run_coverage(_config("mosdepth"))

    assert native.chrom == mos.chrom == _CONTIG
    # Reference length is fixed by the CRAM header -- must be identical.
    assert native.length == mos.length
    # Both count the same CIGAR ops, so aligned bases should match exactly.
    assert native.bases == mos.bases
    # And therefore the reported (2-dp) mean.
    assert native.mean == pytest.approx(mos.mean, abs=0.01)
