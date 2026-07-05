"""
Per-chromosome reductions over the per-base depth vector.

The whole design: `calc_cov(..., per_base=True)` hands us a per-base int32 depth
vector for one chromosome; we immediately reduce it to two compact intermediates
and discard the vector, so peak memory stays ~one chromosome (chr1 ~1 GB
transient) instead of the whole genome.

  depth_histogram()  -> counts per depth value. Sufficient statistic for mean,
                        median, variance, CV, any quantile, breadth. A few KB.
  windowed_means()   -> depth averaged into fixed bins. ~1000x smaller than
                        per-base; feeds every plot and the per-window CN.

Per-base output (d) is handled by `rle_intervals`: coverage is piecewise-constant
between read-block boundaries, so the *change points* of the vector are a
run-length encoding -- orders of magnitude fewer rows than bases.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

# Extreme pileup depths (chrM ~20k, satellite decoys) get clipped into a top bin
# so per-chrom histograms share one length and sum trivially genome-wide. Well
# above any primary-chromosome depth, so the deliverable stats are exact.
DEFAULT_HIST_CAP = 200_000
BREADTH_THRESHOLDS = (1, 5, 10, 15, 20, 30)


@dataclass
class ChromStats:
    n: int                       # positions counted (== chromosome length)
    mean: float
    median: float
    sd: float
    variance: float
    cv: float                    # sd/mean -- coverage uniformity
    mad: float                   # scaled median absolute deviation (robust sd)
    q25: float
    q75: float
    iqr: float
    breadth: dict = field(default_factory=dict)   # depth threshold -> frac positions >= it

    @property
    def robust_cv(self) -> float:
        """MAD/median -- dispersion that ignores the pileup tail sd/mean inflates."""
        return self.mad / self.median if self.median else 0.0


def depth_histogram(base_depth: np.ndarray, cap: int = DEFAULT_HIST_CAP) -> np.ndarray:
    """Counts of positions at each depth 0..cap (depths > cap clipped into cap)."""
    clipped = np.minimum(base_depth, cap)
    return np.bincount(clipped, minlength=cap + 1).astype(np.int64)


def quantile_from_hist(hist: np.ndarray, q: float) -> float:
    """Lower q-quantile of depth, straight from the histogram (no sort, O(cap))."""
    n = int(hist.sum())
    if n == 0:
        return 0.0
    csum = np.cumsum(hist)
    return float(np.searchsorted(csum, q * n, side="left"))


def scaled_mad_from_hist(hist: np.ndarray, median: float, scale: float = 1.4826) -> float:
    """Scaled median absolute deviation: scale * median(|depth - median|).

    Weighted median of the folded deviations, straight from the histogram. The
    1.4826 factor makes it a consistent estimator of sd under normality, so it's
    directly comparable to `sd` but robust to the pileup tail.
    """
    n = int(hist.sum())
    if n == 0:
        return 0.0
    absdev = np.abs(np.arange(hist.size) - median)
    order = np.argsort(absdev, kind="stable")
    csum = np.cumsum(hist[order])
    idx = int(np.searchsorted(csum, 0.5 * n, side="left"))
    return scale * float(absdev[order][idx])


def stats_from_hist(hist: np.ndarray) -> ChromStats:
    depths = np.arange(hist.size, dtype=np.float64)
    n = int(hist.sum())
    if n == 0:
        return ChromStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, {t: 0.0 for t in BREADTH_THRESHOLDS})

    mean = float((depths * hist).sum()) / n
    e_x2 = float((depths * depths * hist).sum()) / n
    variance = max(e_x2 - mean * mean, 0.0)     # clamp tiny negative from fp error
    sd = variance ** 0.5
    median = quantile_from_hist(hist, 0.5)
    mad = scaled_mad_from_hist(hist, median)
    q25 = quantile_from_hist(hist, 0.25)
    q75 = quantile_from_hist(hist, 0.75)

    csum = np.cumsum(hist)
    breadth = {}
    for t in BREADTH_THRESHOLDS:
        ge = n - int(csum[t - 1]) if 0 < t < hist.size else (n if t == 0 else 0)
        breadth[t] = ge / n

    return ChromStats(
        n=n,
        mean=mean,
        median=median,
        sd=sd,
        variance=variance,
        cv=sd / mean if mean else 0.0,
        mad=mad,
        q25=q25,
        q75=q75,
        iqr=q75 - q25,
        breadth=breadth,
    )


def windowed_means(base_depth: np.ndarray, window: int):
    """Mean depth per fixed-size window. Returns (starts, ends, means).

    Uses np.add.reduceat to sum each [start:next_start) segment in one pass; the
    final window is short and divided by its true width.
    """
    length = base_depth.size
    if length == 0:
        return np.array([]), np.array([]), np.array([])
    starts = np.arange(0, length, window)
    sums = np.add.reduceat(base_depth.astype(np.int64), starts)
    ends = np.minimum(starts + window, length)
    means = sums / (ends - starts)
    return starts, ends, means


def rle_intervals(base_depth: np.ndarray, skip_zero: bool = False):
    """Yield (start, end, depth) run-length intervals -> BEDGRAPH rows.

    Change points of the vector = interval boundaries. Interval count scales with
    read-block boundaries, not base count, which is what makes per-base output
    cheap to store (bgzip the result, or convert to BigWig for a browser track).
    """
    if base_depth.size == 0:
        return
    bounds = np.concatenate(([0], np.flatnonzero(np.diff(base_depth)) + 1, [base_depth.size]))
    for i in range(bounds.size - 1):
        start = int(bounds[i])
        depth = int(base_depth[start])
        if skip_zero and depth == 0:
            continue
        yield start, int(bounds[i + 1]), depth


# --- copy number -----------------------------------------------------------

# Autosomes define the diploid (CN=2) reference; exclude sex chroms + mito + alts.
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23)) | frozenset(str(i) for i in range(1, 23))


def is_autosome(chrom: str) -> bool:
    return chrom in AUTOSOMES


def copy_number(depth: float, baseline: float, ploidy: int = 2) -> float:
    """CN ~= ploidy * depth / diploid-baseline. Super approximate: ignores tumor
    purity, ploidy normalization, GC/mappability bias (real callers handle those)."""
    return ploidy * depth / baseline if baseline else 0.0
