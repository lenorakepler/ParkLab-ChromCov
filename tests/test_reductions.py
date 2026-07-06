"""
Unit tests for the reduction/QC classes, on hand-built depth vectors (no CRAM).

These pin the math the pipeline depends on: histogram-based stats (incl. the
lower-quantile convention), windowing, RLE, callability masking (overlap + clip),
and the aneuploidy-aware QC flags.
"""
import numpy as np

from chromcov import qc_flags
from chromcov.depth import ChromDepth, ChromStats, DepthHistogram
from chromcov.qc_flags import copy_number, is_autosome
from chromcov.strata import Strata


def test_histogram_counts_and_stats():
    d = np.array([0, 0, 1, 2, 2, 3], dtype=np.int32)
    cd = ChromDepth(d, cap=10)
    h = cd.histogram().counts
    assert (h[0], h[1], h[2], h[3]) == (2, 1, 2, 1)

    s = cd.stats()
    assert s.n == 6
    assert s.mean == 8 / 6
    # quantile_from_hist uses the lower ("left") quantile convention.
    assert s.median == 1.0
    assert s.q25 == 0.0
    assert s.q75 == 2.0
    assert s.iqr == 2.0
    # breadth = fraction of positions at depth >= threshold.
    assert s.breadth[1] == 4 / 6


def test_depthhistogram_pools_across_chromosomes():
    a = ChromDepth(np.array([2, 2, 2], dtype=np.int32), cap=10).histogram()
    b = ChromDepth(np.array([2, 4], dtype=np.int32), cap=10).histogram()
    pooled = a + b
    assert pooled.n == 5
    assert pooled.counts[2] == 4 and pooled.counts[4] == 1
    assert pooled.stats().mean == (2 + 2 + 2 + 2 + 4) / 5
    # sum() works via __radd__ (0 + hist), for accumulating a genome-wide histogram
    assert isinstance(sum([a, b]), DepthHistogram)
    assert sum([a, b]).n == 5


def test_windowed_means():
    cd = ChromDepth(np.arange(10, dtype=np.int32))
    starts, ends, means = cd.windowed_means(5)
    assert list(starts) == [0, 5]
    assert list(ends) == [5, 10]
    assert means[0] == 2.0 and means[1] == 7.0


def test_windowed_means_ragged_last_window():
    cd = ChromDepth(np.ones(7, dtype=np.int32))
    starts, ends, means = cd.windowed_means(5)
    # last window is 2 bp wide; divided by its true width, still mean 1.0
    assert list(ends) == [5, 7]
    assert list(means) == [1.0, 1.0]


def test_masked_reduces_to_subset():
    cd = ChromDepth(np.array([1, 1, 5, 5], dtype=np.int32), cap=10)
    mask = np.array([True, True, False, False])
    hm = cd.masked(mask).histogram().counts
    assert hm[1] == 2 and hm[5] == 0


def test_strata_mask_overlap_and_clip():
    # intervals [0,2) and [3,10); chrom length 5 -> the second is clipped to 5.
    beds = {"easy": {"chr1": (np.array([0, 3]), np.array([2, 10]))}}
    s = Strata(beds)
    m = s.mask("easy", "chr1", 5)
    assert list(m) == [True, True, False, True, True]
    # a chrom absent from the stratum returns None (not an all-False mask).
    assert s.mask("easy", "chrX", 5) is None
    assert bool(s) and "easy" in s


def test_strata_from_arg_parses_spec(tmp_path):
    bed = tmp_path / "easy.bed"
    bed.write_text("chr1\t0\t2\nchr1\t3\t5\n")
    s = Strata.from_arg(f"easy={bed}")
    assert s.labels() == ["easy"]
    assert list(s.mask("easy", "chr1", 5)) == [True, True, False, True, True]


def test_copy_number_and_autosome():
    assert is_autosome("chr1") and is_autosome("22") and not is_autosome("chrX")
    assert copy_number(20, 10, ploidy=2) == 4.0
    assert copy_number(5, 0) == 0.0  # zero baseline guarded


def _stats(median=10.0, breadth20=0.8, robust_mad=1.0):
    return ChromStats(n=100, mean=median, median=median, sd=1.0, variance=1.0,
                      cv=0.1, mad=robust_mad, q25=median - 1, q75=median + 1, iqr=2.0,
                      breadth={1: 1.0, 10: 0.9, 20: breadth20})


def test_qc_flags_are_aneuploidy_aware():
    st = _stats()
    assert "CN_GAIN" in qc_flags.chrom_flags("chr1", st, cn=3.0, baseline=10)
    assert "CN_LOSS" in qc_flags.chrom_flags("chr1", st, cn=1.0, baseline=10)
    # A single X at CN~1 is normal, not a loss.
    assert qc_flags.chrom_flags("chrX", st, cn=1.0, baseline=10) == []
    # True near-absence (loss of Y) is flagged anywhere.
    assert "CN_DEPLETED" in qc_flags.chrom_flags("chrY", st, cn=0.1, baseline=10)


def test_qc_low_callable_and_uneven():
    assert "LOW_CALLABLE" in qc_flags.chrom_flags("chr1", _stats(breadth20=0.5), cn=2.0, baseline=10)
    assert "UNEVEN" in qc_flags.chrom_flags("chr1", _stats(robust_mad=8.0), cn=2.0, baseline=10)


def test_window_flag_thresholds():
    assert qc_flags.window_flag(3.0) == "GAIN"
    assert qc_flags.window_flag(0.1) == "DEPLETED"
    assert qc_flags.window_flag(1.0) == "LOSS"
    assert qc_flags.window_flag(2.0) == "."


def test_qc_thresholds_are_configurable():
    from chromcov.config import QCThresholds
    st = _stats(median=10.0)
    assert "CN_GAIN" in qc_flags.chrom_flags("chr1", st, cn=2.6, baseline=10)
    # Raise the gain bar above 2.6 -> no longer flagged.
    thr = QCThresholds(gain_cn=3.0)
    assert "CN_GAIN" not in qc_flags.chrom_flags("chr1", st, cn=2.6, baseline=10, thr=thr)
