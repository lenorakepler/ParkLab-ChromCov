"""
On-disk cache of per-contig reduced intermediates (`ReducedContig`), so a re-plot
(or a run that adds contigs) reuses the reduction instead of re-reading and
re-reducing every per-base track -- the reduce pass is the expensive step, the
tracks are already on disk.

One `<chrom>.npz` per contig under a cache dir (default <outdir>/reduced/). Each
file stores a `key` derived from the reduce-relevant config (window size, category
labels, histogram cap, breadth thresholds); a mismatch (or a VERSION bump) makes
`load` return None so the stale entry is recomputed. The cache is disposable:
delete the dir to force a full re-reduce.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np

from .reduce import DepthHistogram
from .result import ReducedContig

VERSION = 1


def cache_key(cfg, categories) -> str:
    """Fingerprint of the inputs that change a ReducedContig (not baseline/ploidy,
    which are applied later at table time)."""
    parts = [str(VERSION), str(cfg.window), str(cfg.hist_cap),
             ",".join(map(str, cfg.breadth_thresholds)),
             ",".join(sorted(categories.labels()))]
    return hashlib.sha1("|".join(parts).encode()).hexdigest()[:16]


def _path(cache_dir, chrom: str) -> Path:
    return Path(cache_dir) / f"{chrom}.npz"


def has(cache_dir, chrom: str) -> bool:
    return _path(cache_dir, chrom).exists()


def cached_chroms(cache_dir) -> list[str]:
    d = Path(cache_dir)
    return sorted(p.stem for p in d.glob("*.npz")) if d.exists() else []


def save(cache_dir, rc: ReducedContig, key: str) -> Path:
    d = Path(cache_dir)
    d.mkdir(parents=True, exist_ok=True)
    labels = list(rc.strata_hist)
    arrays = {
        "key": np.array(key),
        "chrom": np.array(rc.chrom),
        "length": np.array(rc.length),
        "bases": np.array(rc.bases),
        "is_auto": np.array(rc.is_auto),
        "breadth": np.array(rc.hist.breadth_thresholds),
        "hist": np.trim_zeros(rc.hist.counts, "b"),
        "labels": np.array(labels),
        "easy_present": np.array(rc.easy_hist is not None),
        "win_start": np.array([w["start"] for w in rc.win_rows], dtype=np.int64),
        "win_end": np.array([w["end"] for w in rc.win_rows], dtype=np.int64),
        "win_mean": np.array([w["mean"] for w in rc.win_rows], dtype=np.float64),
        "win_easyfrac": np.array([w["easy_frac"] for w in rc.win_rows], dtype=np.float64),
        "win_stratum": np.array([w["stratum"] for w in rc.win_rows]),
    }
    for label in labels:
        arrays[f"sh_{label}"] = np.trim_zeros(rc.strata_hist[label].counts, "b")
        arrays[f"sbp_{label}"] = np.array(rc.strata_bp[label])
    if rc.easy_hist is not None:
        arrays["easy_hist"] = np.trim_zeros(rc.easy_hist.counts, "b")
    dest = _path(d, rc.chrom)
    np.savez_compressed(dest, **arrays)
    return dest


def load(cache_dir, chrom: str, key: str) -> ReducedContig | None:
    """Return the cached ReducedContig, or None if absent or built under a different
    config (stale) -- the caller then re-reduces from the track."""
    path = _path(cache_dir, chrom)
    if not path.exists():
        return None
    with np.load(path, allow_pickle=False) as z:
        if str(z["key"]) != key:
            return None
        breadth = tuple(int(x) for x in z["breadth"])
        hist = DepthHistogram(z["hist"].astype(np.int64), breadth)
        labels = [str(x) for x in z["labels"]]
        strata_hist = {lb: DepthHistogram(z[f"sh_{lb}"].astype(np.int64), breadth) for lb in labels}
        strata_bp = {lb: int(z[f"sbp_{lb}"]) for lb in labels}
        easy_hist = (DepthHistogram(z["easy_hist"].astype(np.int64), breadth)
                     if bool(z["easy_present"]) else None)
        stratum = [str(x) for x in z["win_stratum"]]
        win_rows = [
            {"chrom": chrom, "start": int(s), "end": int(e), "mean": float(m),
             "easy_frac": float(f), "stratum": st}
            for s, e, m, f, st in zip(z["win_start"], z["win_end"], z["win_mean"],
                                      z["win_easyfrac"], stratum)
        ]
        return ReducedContig(
            chrom=chrom, length=int(z["length"]), bases=int(z["bases"]),
            is_auto=bool(z["is_auto"]), hist=hist, strata_hist=strata_hist,
            strata_bp=strata_bp, easy_hist=easy_hist, win_rows=win_rows)
