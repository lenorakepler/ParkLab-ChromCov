"""
Per-chromosome reductions over the per-base depth vector.

Two classes, both here:

`ChromDepth` wraps one chromosome's per-base int32 depth vector (as handed back
by `calc_cov(..., per_base=True)`) and owns every reduction of it:

  .histogram()       -> a DepthHistogram (cached).
  .stats()           -> ChromStats (shortcut for .histogram().stats()).
  .windowed_means()  -> depth averaged into fixed bins. ~1000x smaller than
                        per-base; feeds every plot and the per-window CN.
  .rle_intervals()   -> run-length change-points -> BEDGRAPH rows.
  .masked(mask)      -> a ChromDepth over a boolean-selected subset (callability
                        strata), so a stratum histogram is `.masked(m).histogram()`.

`DepthHistogram` wraps the depth-count array -- a *sufficient statistic* for mean,
median, variance, CV, MAD, any quantile, and breadth-at-depth, in O(max_depth)
with no sort. It adds (`h1 + h2`), so genome-wide and per-stratum stats are just
pooled histograms of the per-chromosome ones -- which is why it isn't tied to a
single ChromDepth.

The design intent: the pipeline reduces each chromosome to compact intermediates
(a few-KB histogram + a windowed track) and discards the vector, so peak memory
stays ~one chromosome (chr1 ~1 GB transient) instead of the whole genome.
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


class DepthHistogram:
    """Counts of positions at each depth 0..cap. A sufficient statistic for the
    whole ChromStats family, and additive: summing per-chromosome histograms
    gives the genome-wide (or per-stratum) one."""

    def __init__(self, counts: np.ndarray, breadth_thresholds=BREADTH_THRESHOLDS):
        self.counts = counts
        self.breadth_thresholds = breadth_thresholds

    @classmethod
    def from_depth(cls, base_depth: np.ndarray, cap: int = DEFAULT_HIST_CAP,
                   breadth_thresholds=BREADTH_THRESHOLDS) -> "DepthHistogram":
        clipped = np.minimum(base_depth, cap)
        counts = np.bincount(clipped, minlength=cap + 1).astype(np.int64)
        return cls(counts, breadth_thresholds)

    @property
    def n(self) -> int:
        return int(self.counts.sum())

    def __add__(self, other) -> "DepthHistogram":
        if other is None or (isinstance(other, int) and other == 0):
            return DepthHistogram(self.counts.copy(), self.breadth_thresholds)
        return DepthHistogram(self.counts + other.counts, self.breadth_thresholds)

    def __radd__(self, other) -> "DepthHistogram":   # enables sum([...])
        return self.__add__(other)

    def quantile(self, q: float) -> float:
        """Lower q-quantile of depth, straight from the histogram (no sort, O(cap))."""
        n = self.n
        if n == 0:
            return 0.0
        csum = np.cumsum(self.counts)
        return float(np.searchsorted(csum, q * n, side="left"))

    def scaled_mad(self, median: float, scale: float = 1.4826) -> float:
        """Scaled median absolute deviation: scale * median(|depth - median|).

        Weighted median of the folded deviations, straight from the histogram. The
        1.4826 factor makes it a consistent estimator of sd under normality, so it's
        directly comparable to `sd` but robust to the pileup tail.
        """
        n = self.n
        if n == 0:
            return 0.0
        absdev = np.abs(np.arange(self.counts.size) - median)
        order = np.argsort(absdev, kind="stable")
        csum = np.cumsum(self.counts[order])
        idx = int(np.searchsorted(csum, 0.5 * n, side="left"))
        return scale * float(absdev[order][idx])

    def stats(self) -> ChromStats:
        counts = self.counts
        depths = np.arange(counts.size, dtype=np.float64)
        n = self.n
        if n == 0:
            return ChromStats(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              {t: 0.0 for t in self.breadth_thresholds})

        mean = float((depths * counts).sum()) / n
        e_x2 = float((depths * depths * counts).sum()) / n
        variance = max(e_x2 - mean * mean, 0.0)     # clamp tiny negative from fp error
        sd = variance ** 0.5
        median = self.quantile(0.5)
        mad = self.scaled_mad(median)
        q25 = self.quantile(0.25)
        q75 = self.quantile(0.75)

        csum = np.cumsum(counts)
        breadth = {}
        for t in self.breadth_thresholds:
            ge = n - int(csum[t - 1]) if 0 < t < counts.size else (n if t == 0 else 0)
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


class ChromDepth:
    """One chromosome's per-base depth vector + every reduction of it."""

    def __init__(self, base_depth: np.ndarray, cap: int = DEFAULT_HIST_CAP,
                 breadth_thresholds=BREADTH_THRESHOLDS):
        self.base_depth = base_depth
        self.cap = cap
        self.breadth_thresholds = breadth_thresholds
        self._hist: DepthHistogram | None = None

    @property
    def length(self) -> int:
        return int(self.base_depth.size)

    def histogram(self) -> DepthHistogram:
        """Depth-count histogram over this chromosome (cached)."""
        if self._hist is None:
            self._hist = DepthHistogram.from_depth(self.base_depth, self.cap, self.breadth_thresholds)
        return self._hist

    def stats(self) -> ChromStats:
        return self.histogram().stats()

    def windowed_means(self, window: int):
        """Mean depth per fixed-size window. Returns (starts, ends, means).

        Uses np.add.reduceat to sum each [start:next_start) segment in one pass;
        the final window is short and divided by its true width.
        """
        length = self.base_depth.size
        if length == 0:
            return np.array([]), np.array([]), np.array([])
        starts = np.arange(0, length, window)
        sums = np.add.reduceat(self.base_depth.astype(np.int64), starts)
        ends = np.minimum(starts + window, length)
        means = sums / (ends - starts)
        return starts, ends, means

    def rle_intervals(self, skip_zero: bool = False):
        """Yield (start, end, depth) run-length intervals -> BEDGRAPH rows.

        Change points of the vector = interval boundaries. Interval count scales
        with read-block boundaries, not base count, which is what makes per-base
        output cheap to store (bgzip the result, or convert to BigWig).
        """
        d = self.base_depth
        if d.size == 0:
            return
        bounds = np.concatenate(([0], np.flatnonzero(np.diff(d)) + 1, [d.size]))
        for i in range(bounds.size - 1):
            start = int(bounds[i])
            depth = int(d[start])
            if skip_zero and depth == 0:
                continue
            yield start, int(bounds[i + 1]), depth

    def masked(self, mask: np.ndarray) -> "ChromDepth":
        """A ChromDepth over the boolean-selected subset (e.g. a callability tier)."""
        return ChromDepth(self.base_depth[mask], self.cap, self.breadth_thresholds)


# --- copy number -----------------------------------------------------------

# Autosomes define the diploid (CN=2) reference; exclude sex chroms + mito + alts.
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23)) | frozenset(str(i) for i in range(1, 23))


def is_autosome(chrom: str) -> bool:
    return chrom in AUTOSOMES


def copy_number(depth: float, baseline: float, ploidy: int = 2) -> float:
    """CN ~= ploidy * depth / diploid-baseline. Super approximate: ignores tumor
    purity, ploidy normalization, GC/mappability bias (real callers handle those)."""
    return ploidy * depth / baseline if baseline else 0.0
