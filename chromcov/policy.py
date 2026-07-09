"""
Coverage-QC policy: what counts as normal.

Consolidates the scattered "normal" decisions (was qc_flags.py +
QCReport.baseline/finalize):

  - the STANDARD (autosomal) set that defines the diploid reference,
  - the diploid `baseline` choice,
  - approximate `copy_number` (the ratio metric),
  - the per-chromosome / per-window abnormality flags,
  - `finalize`, which folds baseline + copy number + flags into the RunResult's
    per-chromosome rows.

These are coverage-QC *abnormality* flags (CN_GAIN, LOW_DEPTH, ...), not SAM read
flags -- those, and the read filtering they drive, live in `chromcov.filtering`.
Everything here operates on reduced stats (a ChromStats or a scalar), never on a
per-base vector, so there are no numpy-heavy reductions -- just thresholds.
"""
from __future__ import annotations

from .config.schema import QCThresholds
from .present.frames import ChromCoverage
from .reduce import ChromStats

# Autosomes define the diploid (CN=2) reference; exclude sex chroms + mito + alts.
# This is the STANDARD contig set the baseline is computed over.
AUTOSOMES = frozenset(f"chr{i}" for i in range(1, 23)) | frozenset(str(i) for i in range(1, 23))


def is_autosome(chrom: str) -> bool:
    return chrom in AUTOSOMES


# Alias with the plan's vocabulary: the "standard set" membership test.
def is_standard(chrom: str) -> bool:
    return is_autosome(chrom)


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


def baseline(result) -> tuple[float, str]:
    """Diploid (CN=2) reference depth for a RunResult (cached on it). Prefer the
    callability-masked easy-autosomal median (drops segdup/centromere pileups);
    fall back to all autosomes, then to a length-weighted mean if the subset has
    no autosomes."""
    if result._baseline is not None:
        return result._baseline, result._baseline_source

    val, src = 0.0, ""
    if result.cfg.baseline == "easy-autosomal-median" and result.easy_autosomal_hist is not None:
        val, src = result.easy_autosomal_hist.quantile(0.5), "easy-autosomal median"
    if not val and result.autosomal_hist is not None:
        val, src = result.autosomal_hist.quantile(0.5), "autosomal median"
    if not val:
        total_len = sum(result.lengths.values()) or 1
        val = sum(result.per_chrom_stats[c].mean * result.lengths[c]
                  for c in result.chroms) / total_len or 1.0
        src = "length-weighted mean (no autosomes in subset)"

    result._baseline, result._baseline_source = val, src
    return val, src


def finalize(result) -> None:
    """Fold baseline + copy number + QC flags into the RunResult's per-chromosome
    rows (and its flagged list). Pure over the accumulated result -- no I/O."""
    base, _ = baseline(result)
    ploidy = result.cfg.ploidy

    result.rows = []
    result.flagged = []
    for chrom in result.chroms:
        s = result.per_chrom_stats[chrom]
        cn = copy_number(s.mean, base, ploidy)
        fl = chrom_flags(chrom, s, cn, base, result.cfg.qc)
        if fl:
            result.flagged.append((chrom, fl))
        result.rows.append(ChromCoverage(
            chrom=chrom, length=result.lengths[chrom], bases=result.bases[chrom],
            stats=s, copy_number=cn, flags=fl))
