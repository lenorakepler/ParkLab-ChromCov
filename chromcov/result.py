"""
RunResult -- the accumulator a pipeline run folds into (replaces QCReport's ~15
mutable attributes).

It holds the raw accumulated intermediates (per-chrom stats, pooled + per-category
histograms, windowed rows) and knows how to reduce one contig's per-base vector
into them (`reduce_chrom`) or record one contig's mean total (`add_mean`). The
policy tail (`policy.baseline` / `policy.finalize`) and the presentation tail
(`present.frames`) are pure functions over this object -- it does no I/O itself.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .categories import Strata
from .policy import is_autosome
from .present.frames import ChromCoverage
from .reduce import ChromDepth, ChromStats, DepthHistogram


@dataclass
class RunResult:
    cfg: object
    depth: object = None                     # a pipeline.Depth (kept untyped to avoid a cycle)
    categories: Strata | None = None
    preflight: dict | None = None

    chroms: list[str] = field(default_factory=list)
    lengths: dict[str, int] = field(default_factory=dict)
    bases: dict[str, int] = field(default_factory=dict)
    per_chrom_stats: dict[str, ChromStats] = field(default_factory=dict)
    win_rows: list[dict] = field(default_factory=list)

    autosomal_hist: DepthHistogram | None = None
    easy_autosomal_hist: DepthHistogram | None = None
    strata_hist: dict[str, DepthHistogram] = field(default_factory=dict)
    strata_bp: dict[str, int] = field(default_factory=dict)

    # filled by policy.finalize
    rows: list[ChromCoverage] = field(default_factory=list)
    flagged: list[tuple[str, list[str]]] = field(default_factory=list)
    _baseline: float | None = None
    _baseline_source: str = ""

    def __post_init__(self):
        if self.categories is None:
            self.categories = Strata({})
        if not self.strata_bp:
            self.strata_bp = {label: 0 for label in self.categories.labels()}

    # --- accumulation (the "middle") ---------------------------------------

    def add_mean(self, chrom: str, length: int, bases: int) -> None:
        """Record one contig's mean total (the track-free MEAN path)."""
        self.chroms.append(chrom)
        self.lengths[chrom] = length
        self.bases[chrom] = bases

    def reduce_chrom(self, chrom: str, length: int, base_depth: np.ndarray) -> None:
        """
        Reduce one contig's per-base vector to hist + windows (+ categories),
        accumulate the pooled histograms, then drop the vector (the FULL path).
        """
        total_depth = int(base_depth.sum())
        depth = ChromDepth(base_depth, cap=self.cfg.hist_cap,
                           breadth_thresholds=self.cfg.breadth_thresholds)

        hist = depth.histogram()
        self.chroms.append(chrom)
        self.per_chrom_stats[chrom] = hist.stats()
        self.lengths[chrom] = length
        self.bases[chrom] = total_depth

        if is_autosome(chrom):
            self.autosomal_hist = hist if self.autosomal_hist is None else self.autosomal_hist + hist

        masks: dict[str, np.ndarray] = {}
        for label in self.categories.labels():
            mask = self.categories.mask(label, chrom, length)
            if mask is None:
                continue
            masks[label] = mask
            h = depth.masked(mask).histogram()
            self.strata_hist[label] = h if label not in self.strata_hist else self.strata_hist[label] + h
            self.strata_bp[label] += int(mask.sum())
            if label == "easy" and is_autosome(chrom):
                self.easy_autosomal_hist = h if self.easy_autosomal_hist is None else self.easy_autosomal_hist + h

        starts, ends, means = depth.windowed_means(self.cfg.window)
        widths = ends - starts
        # Per-window fraction in each category, plus the dominant tier -- the scatter
        # colors every window by its tier (easy/difficult/extreme) instead of only
        # keeping callable ones. `easy_frac` is retained for the optional
        # callable-only view; with no categories a window is easy_frac=1, stratum="".
        fracs = {label: np.add.reduceat(m.astype(np.int64), starts) / widths
                 for label, m in masks.items()}
        easy_ef = fracs.get("easy", np.ones_like(means))
        if fracs:
            order = list(fracs)
            dom = [order[i] for i in np.vstack([fracs[label] for label in order]).argmax(axis=0).tolist()]
        else:
            dom = [""] * len(means)
        for s, e, m, f, st in zip(starts.tolist(), ends.tolist(), means.tolist(),
                                  easy_ef.tolist(), dom):
            self.win_rows.append({"chrom": chrom, "start": s, "end": e, "mean": m,
                                  "easy_frac": f, "stratum": st})

    # --- views over the accumulated state ----------------------------------

    @property
    def baseline(self) -> tuple[float | None, str]:
        """The cached diploid baseline (value, source); computed by policy.finalize."""
        return self._baseline, self._baseline_source

    def coverage_rows(self) -> list[ChromCoverage]:
        """The per-chromosome rows: the finalized rows (with stats/CN/flags) for a
        --full run, else plain mean-only rows for the MEAN path."""
        if self.rows:
            return self.rows
        return [ChromCoverage(chrom=c, length=self.lengths[c], bases=self.bases[c])
                for c in self.chroms]
