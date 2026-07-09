"""
The reporting layer (polars): the per-chromosome frame carries only the columns
the rows populate (base for mean-only; + stats/CN/flags for --full), and the
windowed frame derives copy number + focal flag vectorized.
"""
import numpy as np

from chromcov.reduce import ChromDepth
from chromcov.present.frames import (
    ChromCoverage,
    coverage_frame,
    windows_frame,
    write_table,
)


def test_mean_only_frame_is_base_columns():
    rows = [ChromCoverage(chrom="chr1", length=100, bases=1500)]
    assert rows[0].mean == 15.0
    df = coverage_frame(rows)
    assert df.columns == ["chrom", "length", "bases", "mean_coverage"]
    assert df["mean_coverage"][0] == 15.0


def test_full_frame_carries_stats_and_cn():
    stats = ChromDepth(np.array([10, 10, 20, 20], dtype=np.int32), cap=50).stats()
    rows = [ChromCoverage(chrom="chr1", length=4, bases=60, stats=stats,
                          copy_number=1.5, flags=["CN_LOSS"])]
    df = coverage_frame(rows)
    assert {"median", "copy_number", "flags"} <= set(df.columns)
    assert df["copy_number"][0] == 1.5
    assert df["flags"][0] == "CN_LOSS"


def test_empty_flags_render_ok():
    stats = ChromDepth(np.array([5, 5, 5], dtype=np.int32), cap=50).stats()
    rows = [ChromCoverage(chrom="chr1", length=3, bases=15, stats=stats,
                          copy_number=2.0, flags=[])]
    assert coverage_frame(rows)["flags"][0] == "OK"


def test_windows_frame_vectorized_cn_and_flag():
    win = [
        {"chrom": "chr1", "start": 0, "end": 10, "mean": 20.0, "easy_frac": 1.0},   # cn 4 -> GAIN
        {"chrom": "chr1", "start": 10, "end": 20, "mean": 5.0, "easy_frac": 1.0},   # cn 1 -> LOSS
        {"chrom": "chr1", "start": 20, "end": 30, "mean": 10.0, "easy_frac": 1.0},  # cn 2 -> .
    ]
    df = windows_frame(win, baseline=10.0, ploidy=2)
    assert df.columns == ["chrom", "start", "end", "mean_depth", "copy_number", "flag"]
    assert df["copy_number"].to_list() == [4.0, 1.0, 2.0]
    assert df["flag"].to_list() == ["GAIN", "LOSS", "."]


def test_write_table_tsv_header_and_values(tmp_path):
    rows = [ChromCoverage(chrom="chr1", length=100, bases=1500)]
    path = tmp_path / "coverage.tsv"
    write_table(coverage_frame(rows), path)
    lines = path.read_text().splitlines()
    assert lines[0].split("\t") == ["chrom", "length", "bases", "mean_coverage"]
    assert lines[1].split("\t")[0] == "chr1"
