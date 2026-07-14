"""
The per-contig reduce cache: a ReducedContig round-trips through disk unchanged,
folding a cached contig gives the same pooled result as folding the fresh one, and
the cache key invalidates when the reduce-relevant config changes.
"""
from types import SimpleNamespace

import numpy as np

from chromcov import reduce_cache
from chromcov.categories import Strata
from chromcov.result import RunResult


def _cfg(window=5):
    return SimpleNamespace(window=window, hist_cap=1000, breadth_thresholds=(1, 10, 20))


class _StubStrata:
    """Categories stub: one 'easy' tier masking the first half of the contig."""
    def labels(self):
        return ["easy"]

    def mask(self, label, chrom, length):
        m = np.zeros(length, dtype=bool)
        m[: length // 2] = True
        return m


def _depth():
    # chr1 so it counts as an autosome (drives the easy-autosomal baseline path)
    return np.array([10, 10, 12, 8, 30, 0, 5, 5, 20, 20], dtype=np.int32)


def test_roundtrip_matches_fresh(tmp_path):
    res = RunResult(cfg=_cfg(), categories=_StubStrata())
    rc = res._reduce_one("chr1", 10, _depth())

    key = reduce_cache.cache_key(res.cfg, res.categories)
    reduce_cache.save(tmp_path, rc, key)
    loaded = reduce_cache.load(tmp_path, "chr1", key)

    assert loaded.chrom == "chr1" and loaded.length == 10 and loaded.bases == rc.bases
    assert loaded.is_auto is True
    assert loaded.hist.stats().mean == rc.hist.stats().mean
    assert loaded.win_rows == rc.win_rows
    assert set(loaded.strata_hist) == {"easy"}
    assert loaded.strata_bp == rc.strata_bp
    assert loaded.easy_hist is not None


def test_fold_of_cached_equals_fold_of_fresh(tmp_path):
    fresh = RunResult(cfg=_cfg(), categories=_StubStrata())
    rc = fresh._reduce_one("chr1", 10, _depth())
    fresh.fold(rc)

    key = reduce_cache.cache_key(fresh.cfg, fresh.categories)
    reduce_cache.save(tmp_path, rc, key)
    cached = RunResult(cfg=_cfg(), categories=_StubStrata())
    cached.fold(reduce_cache.load(tmp_path, "chr1", key))

    assert cached.per_chrom_stats["chr1"].mean == fresh.per_chrom_stats["chr1"].mean
    assert cached.win_rows == fresh.win_rows
    assert cached.strata_bp == fresh.strata_bp
    assert cached.autosomal_hist.quantile(0.5) == fresh.autosomal_hist.quantile(0.5)
    assert cached.easy_autosomal_hist.quantile(0.5) == fresh.easy_autosomal_hist.quantile(0.5)


def test_key_invalidates_on_config_change(tmp_path):
    res = RunResult(cfg=_cfg(window=5), categories=_StubStrata())
    rc = res._reduce_one("chr1", 10, _depth())
    key = reduce_cache.cache_key(res.cfg, res.categories)
    reduce_cache.save(tmp_path, rc, key)

    other_key = reduce_cache.cache_key(_cfg(window=7), _StubStrata())
    assert other_key != key
    assert reduce_cache.load(tmp_path, "chr1", other_key) is None      # stale -> miss
    assert reduce_cache.load(tmp_path, "chr1", key) is not None        # matching -> hit


def test_cached_chroms_lists_saved(tmp_path):
    res = RunResult(cfg=_cfg(), categories=Strata({}))
    key = reduce_cache.cache_key(res.cfg, res.categories)
    for c in ("chr1", "chr2"):
        reduce_cache.save(tmp_path, res._reduce_one(c, 10, _depth()), key)
    assert reduce_cache.cached_chroms(tmp_path) == ["chr1", "chr2"]
