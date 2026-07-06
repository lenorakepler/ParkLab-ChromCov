"""
Copy number + abnormality flags for coverage QC.

These are coverage-QC *abnormality* flags -- human-readable strings like CN_GAIN,
CN_LOSS, LOW_DEPTH, UNEVEN, applied to a chromosome or a window. They are NOT SAM
read flags (unmapped/duplicate/...); those, and the read filtering they drive,
live in `config.py` (`SAM_FLAGS`, `ReadFilter`). The module is named `qc_flags`
to keep the two senses of "flag" distinct.

Approximate copy number (`CN ~= ploidy * depth / diploid-baseline`) and the
per-chromosome / per-window flags derived from it. Kept together because both
operate on reduced stats (a ChromStats or a scalar depth), never on the per-base
vector -- so this module has no numpy-heavy reductions, just thresholds.

The copy-number logic is aneuploidy-aware: gains/losses are only auto-flagged on
autosomes (a male X at CN~1 is normal, not a loss), while true near-absence
(CN~0, e.g. loss of Y) is flagged anywhere. Dispersion uses the robust MAD/median
(`stats.robust_cv`), not sd/mean, so the multi-mapping pileup tail doesn't trip a
false UNEVEN flag.
"""
from __future__ import annotations

from .config import QCThresholds
from .depth import ChromStats

# Autosomes define the diploid (CN=2) reference; exclude sex chroms + mito + alts.
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23)) | frozenset(str(i) for i in range(1, 23))

def is_autosome(chrom: str) -> bool:
    return chrom in AUTOSOMES

def copy_number(depth: float, baseline: float, ploidy: int = 2) -> float:
    """CN ~= ploidy * depth / diploid-baseline. Super approximate: ignores tumor
    purity, ploidy normalization, GC/mappability bias (real callers handle those)."""
    return ploidy * depth / baseline if baseline else 0.0

def chrom_flags(chrom: str, stats: ChromStats, cn: float,
                baseline: float, thr: QCThresholds = QCThresholds()) -> list[str]:
    flags: list[str] = []
    autosome = is_autosome(chrom)

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
