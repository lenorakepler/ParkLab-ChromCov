"""
The combined coverage table: one ChromCoverage row type, dynamic columns.
--fast rows carry only mean; full rows carry stats (+ copy number / flags).
"""
import io

import numpy as np

from chromcov.analysis import ChromDepth
from chromcov.result import BASE_COLUMNS, columns, write_rows
from chromcov.result import ChromCoverage


def test_fast_row_is_base_columns_only():
    rows = [ChromCoverage(chrom="chr1", length=100, bases=1500)]
    assert columns(rows) == BASE_COLUMNS
    r = rows[0].as_row()
    assert r["mean_coverage"] == 15.0 and "median" not in r


def test_full_row_carries_stats_and_cn():
    stats = ChromDepth(np.array([10, 10, 20, 20], dtype=np.int32), cap=50).stats()
    rows = [ChromCoverage(chrom="chr1", length=4, bases=60, stats=stats,
                          copy_number=1.5, flags=["CN_LOSS"])]
    cols = columns(rows)
    assert "median" in cols and "copy_number" in cols and "flags" in cols
    r = rows[0].as_row()
    assert r["copy_number"] == 1.5 and r["flags"] == "CN_LOSS"


def test_write_rows_header_matches_columns():
    stats = ChromDepth(np.array([5, 5, 5], dtype=np.int32), cap=50).stats()
    rows = [ChromCoverage(chrom="chr1", length=3, bases=15, stats=stats, copy_number=2.0, flags=[])]
    buf = io.StringIO()
    write_rows(rows, buf)
    header = buf.getvalue().splitlines()[0].split("\t")
    assert header == columns(rows)
    assert "OK" in buf.getvalue()   # empty flags -> "OK"
