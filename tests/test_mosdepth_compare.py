"""
Cross-check the native calculator against the optional mosdepth add-on.

mosdepth is no longer an in-tool backend -- it's a standalone script
(scripts/mosdepth_coverage.py) that emits the same coverage.tsv format. This
test runs both on one contig and asserts they agree, which is exactly the
validation the old two-backend dispatch existed for, minus the machinery.

Both count aligned bases over M/=/X CIGAR ops with no overlap correction and
`mean = bases / length`, and the add-on defaults to chromcov's pinned flag mask,
so the numbers should match to rounding.

    uv run --extra dev pytest tests/test_mosdepth_compare.py -v

Skips cleanly when mosdepth or the data/ CRAM aren't present.
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from chromcov import Depth, run
from chromcov.config.schema import Config

_REPO_ROOT = Path(__file__).resolve().parents[1]
_DATA = _REPO_ROOT / "data"
_CRAM = _DATA / "COLO829T_TEST.cram"
_REF = _DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
_SCRIPT = _REPO_ROOT / "scripts" / "mosdepth_coverage.py"

_CONTIG = "chr21"  # smallest standard autosome; CNV-interesting in COLO829T

requires_mosdepth = pytest.mark.skipif(
    shutil.which("mosdepth") is None, reason="mosdepth not on PATH"
)
requires_data = pytest.mark.skipif(
    not (_CRAM.exists() and _REF.exists()), reason="data/ CRAM or reference missing"
)


@requires_data
@requires_mosdepth
def test_native_matches_mosdepth_addon(tmp_path):
    (native,) = run(Config(cram=_CRAM, reference=_REF, chroms=(_CONTIG,)),
                    depth=Depth.MEAN).coverage_rows()

    subprocess.run(
        [sys.executable, str(_SCRIPT), "--cram", str(_CRAM), "--reference", str(_REF),
         "--chrom", _CONTIG, "--outdir", str(tmp_path), "--no-plots"],
        check=True,
    )
    with (tmp_path / "coverage.tsv").open() as fh:
        rows = {r["chrom"]: r for r in csv.DictReader(fh, delimiter="\t")}
    mos = rows[_CONTIG]

    assert int(mos["length"]) == native.length     # fixed by the CRAM header
    assert int(mos["bases"]) == native.bases        # same CIGAR-op accounting
    assert float(mos["mean_coverage"]) == pytest.approx(native.mean, abs=0.01)
