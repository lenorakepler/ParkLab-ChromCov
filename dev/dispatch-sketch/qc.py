"""
Abnormality flags for coverage QC.

Turns the per-chromosome stats into a short list of human-readable flags, so a
reviewer sees *what's off* without reading the whole table. Thresholds are
configurable (QCThresholds) and the copy-number logic is aneuploidy-aware:
gains/losses are only auto-flagged on autosomes (a male X at CN~1 is normal, not
a loss), while true near-absence (CN~0, e.g. loss of Y) is flagged anywhere.

Dispersion uses the robust MAD/median (`stats.robust_cv`), not sd/mean, so the
multi-mapping pileup tail doesn't trip a false UNEVEN flag.
"""
from __future__ import annotations

from dataclasses import dataclass

import analysis


@dataclass
class QCThresholds:
    gain_cn: float = 2.5          # >= -> CN_GAIN (autosomes)
    loss_cn: float = 1.5          # <= -> CN_LOSS (autosomes)
    depleted_cn: float = 0.25     # <= -> CN_DEPLETED (any contig; e.g. LOY, homozygous del)
    min_median: float = 10.0      # median depth below this -> LOW_DEPTH
    uneven_robust_cv: float = 0.5  # MAD/median above this -> UNEVEN (focal CNV / artifact)
    min_breadth_20x: float = 0.70  # callable breadth below this -> LOW_CALLABLE
    extreme_median_mult: float = 5.0  # median > mult * baseline -> EXTREME_DEPTH (artifact/mito)


def chrom_flags(chrom: str, stats: analysis.ChromStats, cn: float,
                baseline: float, thr: QCThresholds = QCThresholds()) -> list[str]:
    flags: list[str] = []
    autosome = analysis.is_autosome(chrom)

    # Copy number
    if cn <= thr.depleted_cn:
        flags.append("CN_DEPLETED")
    elif autosome and cn >= thr.gain_cn:
        flags.append("CN_GAIN")
    elif autosome and cn <= thr.loss_cn:
        flags.append("CN_LOSS")

    # Depth / callability
    if stats.median < thr.min_median:
        flags.append("LOW_DEPTH")
    if baseline and stats.median > thr.extreme_median_mult * baseline:
        flags.append("EXTREME_DEPTH")
    if stats.breadth.get(20, 1.0) < thr.min_breadth_20x:
        flags.append("LOW_CALLABLE")

    # Uniformity (robust)
    if stats.robust_cv > thr.uneven_robust_cv:
        flags.append("UNEVEN")

    return flags


def window_flag(cn: float, thr: QCThresholds = QCThresholds()) -> str:
    """Focal per-window call, for spotting intrachromosomal segments/breakpoints."""
    if cn <= thr.depleted_cn:
        return "DEPLETED"
    if cn >= thr.gain_cn:
        return "GAIN"
    if cn <= thr.loss_cn:
        return "LOSS"
    return "."
