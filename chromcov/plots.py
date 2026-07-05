"""
Plots, read from the saved windowed track (never from a live CRAM recompute) so
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

# Primary assembly in karyotypic order; decoys/unplaced omitted from the headline
# plots (their per-base means are multi-mapping artifacts).
PRIMARY = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


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
    """windows: rows with chrom/start/mean(/easy_frac), restricted to PRIMARY. Y is
    copy ratio (2*mean/baseline), capped so pileup spikes don't flatten the axis.

    min_easy_frac > 0 drops windows that are mostly outside the callable ('easy')
    mask -- i.e. the centromere/segdup pileups -- so real segmental CN steps stand
    out instead of vertical repeat spikes.
    """
    fig, ax = plt.subplots(figsize=(16, 4))

    offset = 0
    xticks, xlabels = [], []
    for chrom in PRIMARY:
        rows = [w for w in windows
                if w["chrom"] == chrom and w.get("easy_frac", 1.0) >= min_easy_frac]
        if not rows:
            continue
        xs = [offset + w["start"] for w in rows]
        ys = [min(ploidy * w["mean"] / baseline, cap_cn) if baseline else 0 for w in rows]
        ax.scatter(xs, ys, s=2, alpha=0.4)
        span = max(w["start"] for w in rows)
        xticks.append(offset + span / 2)
        xlabels.append(chrom.replace("chr", ""))
        offset += span + 1

    ax.axhline(ploidy, color="grey", ls="--", lw=1)   # CN=2 diploid line
    ax.set_ylim(0, cap_cn)
    ax.set_ylabel(f"approx copy number (cap {cap_cn:g})")
    masked = f", callable windows only (easy≥{min_easy_frac:g})" if min_easy_frac else ""
    ax.set_title(f"Windowed copy number across the genome{masked}")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=120)
    plt.close(fig)
    return out_path
