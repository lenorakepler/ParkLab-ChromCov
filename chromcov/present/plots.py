"""
Plots, read from the saved windowed means (never from a live CRAM recompute) so
styling is decoupled from the expensive coverage pass.

  bar_by_chromosome  -> mean depth per primary chromosome (the headline table).
  scatter_windows    -> windowed copy-ratio along the genome; intrachromosomal
                        CNV breakpoints show up as step changes between windows.

matplotlib is imported with the Agg backend so this runs headless (CI, servers).
"""
from __future__ import annotations

from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ..categories import STRATUM_ORDER

# Primary assembly in karyotypic order; decoys/unplaced omitted from the headline
# plots (their per-base means are multi-mapping artifacts).
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

# Callability tier -> point color (best -> worst callability), + the fallback
# color when no strata were supplied.
STRATA_COLORS = {"easy": "#2ca02c", "difficult": "#ff7f0e", "extreme": "#d62728"}
UNSTRATIFIED_COLOR = "#4C72B0"


def bar_by_chromosome(chrom_means: dict[str, float], out_path: Path, baseline: float | None = None):
    chroms = [c for c in PRIMARY if c in chrom_means]
    vals = [chrom_means[c] for c in chroms]

    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(chroms, vals, color="#4C72B0")
    if baseline:
        ax.axhline(baseline, color="crimson", ls="--", lw=1, label=f"autosomal median ({baseline:.1f}x)")
        ax.legend()
    ax.set_ylabel("mean depth (x)")
    ax.set_title("Mean coverage per chromosome")
    ax.tick_params(axis="x", rotation=90)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path


def scatter_windows(windows: list[dict], out_path: Path, baseline: float, ploidy: int = 2,
                    cap_cn: float = 6.0, min_easy_frac: float = 0.0):
    """windows: rows with chrom/start/mean(/easy_frac/stratum), restricted to
    PRIMARY. Y is copy ratio (ploidy*mean/baseline), capped so pileup spikes don't
    flatten the axis.

    When strata were supplied, every window is shown and colored by its dominant
    callability tier (easy/difficult/extreme), so segmental CN steps and the
    repeat/centromere pileups are visible *and* distinguishable rather than the
    latter being silently dropped. `min_easy_frac > 0` instead restricts to the
    callable ('easy') windows (the old callable-only view).
    """
    fig, ax = plt.subplots(figsize=(16, 4))

    offset = 0
    xticks, xlabels, boundaries = [], [], []
    by_tier: dict[str, tuple[list, list]] = {}   # tier -> (xs, ys)
    for chrom in PRIMARY:
        rows = [w for w in windows
                if w["chrom"] == chrom and w.get("easy_frac", 1.0) >= min_easy_frac]
        if not rows:
            continue
        boundaries.append(offset)   # left edge of this chromosome's band
        for w in rows:
            x = offset + w["start"]
            y = min(ploidy * w["mean"] / baseline, cap_cn) if baseline else 0
            tier = w.get("stratum", "") or ""
            by_tier.setdefault(tier, ([], []))
            by_tier[tier][0].append(x)
            by_tier[tier][1].append(y)
        span = max(w["start"] for w in rows)
        xticks.append(offset + span / 2)
        xlabels.append(chrom.replace("chr", ""))
        offset += span + 1

    # Vertical dividers between chromosome bands (behind the points).
    for b in boundaries[1:]:
        ax.axvline(b, color="0.8", lw=0.5, zorder=0)

    stratified = any(tier for tier in by_tier)
    if stratified:
        order = [t for t in STRATUM_ORDER if t in by_tier] + \
                [t for t in by_tier if t and t not in STRATUM_ORDER]
        for tier in order:
            xs, ys = by_tier[tier]
            ax.scatter(xs, ys, s=2, alpha=0.4, color=STRATA_COLORS.get(tier, "#888888"), label=tier)
        ax.legend(title="callability", markerscale=4, fontsize=7, loc="upper right")
    else:
        xs, ys = by_tier.get("", ([], []))
        ax.scatter(xs, ys, s=2, alpha=0.4, color=UNSTRATIFIED_COLOR)

    ax.axhline(ploidy, color="grey", ls="--", lw=1)   # CN=2 diploid line
    ax.set_ylim(0, cap_cn)
    ax.set_ylabel(f"approx copy number (cap {cap_cn:g})")
    if min_easy_frac:
        note = f" — callable windows only (easy≥{min_easy_frac:g})"
    elif stratified:
        note = " — colored by callability tier"
    else:
        note = ""
    ax.set_title(f"Windowed copy number across the genome{note}")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path
