"""
Coverage analysis driver: per-chromosome stats + windowed track + copy number +
plots, all from one memory-bounded pass over the native per-base depth.

    uv run --with matplotlib python dev/dispatch-sketch/analyze.py --chroms chr20,chr21,chrX,chrY,chrM
    uv run --with matplotlib python dev/dispatch-sketch/analyze.py            # whole genome

Outputs (under out/ by default):
  coverage.stats.tsv    per-chromosome mean/median/sd/CV/IQR/breadth/CN
  coverage.windows.bed  windowed depth + copy number (feeds the plots)
  coverage.perbase.bedgraph[.gz]   optional (--per-base), RLE per-base depth
  coverage.bar.png / coverage.scatter.png

Runs the native backend directly (per_base=True) so this needs no mosdepth.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT))

from chromcov.read_filter import ReadFilter, calc_cov
import analysis
import plots
import qc
import strata as strata_mod
import validate
from config import CoverageConfig

_DATA = _REPO_ROOT / "data"
_CRAM = _DATA / "COLO829T_TEST.cram"
_REF = _DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

STATS_COLUMNS = [
    "chrom", "length", "mean", "median", "sd", "cv", "mad", "robust_cv",
    "q25", "q75", "iqr", "breadth_1x", "breadth_10x", "breadth_20x",
    "copy_number", "flags",
]


def _write_stats(rows: list[dict], path: Path):
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STATS_COLUMNS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _write_windows(windows: list[dict], path: Path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["chrom", "start", "end", "mean_depth", "copy_number", "flag"])
        for x in windows:
            w.writerow([x["chrom"], x["start"], x["end"], f"{x['mean']:.2f}",
                        f"{x['cn']:.3f}", x.get("flag", ".")])


STRATA_COLUMNS = ["stratum", "region_bp", "pct_of_analyzed", "mean", "median", "sd", "cv",
                  "breadth_1x", "breadth_10x", "breadth_20x"]


def _write_strata(strata_hist: dict, strata_bp: dict, path: Path):
    total_bp = sum(strata_bp.values()) or 1
    order = [s for s in strata_mod.STRATUM_ORDER if s in strata_hist] + \
            [s for s in strata_hist if s not in strata_mod.STRATUM_ORDER]
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=STRATA_COLUMNS, delimiter="\t")
        w.writeheader()
        for label in order:
            s = analysis.stats_from_hist(strata_hist[label])
            w.writerow({
                "stratum": label, "region_bp": strata_bp[label],
                "pct_of_analyzed": round(100 * strata_bp[label] / total_bp, 2),
                "mean": round(s.mean, 2), "median": round(s.median, 2),
                "sd": round(s.sd, 2), "cv": round(s.cv, 3),
                "breadth_1x": round(s.breadth[1], 4), "breadth_10x": round(s.breadth[10], 4),
                "breadth_20x": round(s.breadth[20], 4),
            })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default="", help="comma-separated subset (default: all)")
    ap.add_argument("--window", type=int, default=10_000)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--per-base", action="store_true", help="also write RLE bedgraph (bgzip if .gz)")
    ap.add_argument("--strata", default="",
                    help="callability strata as label=bed[,label=bed]; e.g. "
                         "easy=SMaHT_easy_hg38.bed.gz,difficult=...,extreme=...")
    ap.add_argument("--outdir", default=str(_REPO_ROOT / "out"))
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Fail fast on bad inputs before spending minutes decoding the CRAM.
    cfg = CoverageConfig(cram=_CRAM, reference=_REF, min_mapping_quality=args.min_mapq)
    report = validate.preflight(cfg)
    print(f"[preflight] ok: sorted, indexed, reference {report['reference_check']['status']}")

    rf = ReadFilter(min_mapping_quality=args.min_mapq)
    cram = pysam.AlignmentFile(str(_CRAM), "rc", reference_filename=str(_REF),
                               index_filename=str(_CRAM) + ".crai")
    chroms = args.chroms.split(",") if args.chroms else list(cram.references)

    # First pass: reduce each chromosome to hist + windows, accumulate the
    # autosomal pooled histogram for the CN baseline. Per-base vector discarded
    # immediately -> peak memory ~one chromosome.
    strata_beds = {label: strata_mod.load_bed(path)
                   for label, path in strata_mod.parse_strata_arg(args.strata).items()}

    per_chrom_stats: dict[str, analysis.ChromStats] = {}
    lengths: dict[str, int] = {}
    win_rows: list[dict] = []
    autosomal_hist = None
    # Pooled genome-wide histogram per stratum, plus the easy-autosomal one used
    # as the callability-masked diploid baseline.
    strata_hist: dict[str, np.ndarray] = {}
    strata_bp: dict[str, int] = {label: 0 for label in strata_beds}
    easy_autosomal_hist = None
    perbase_fh = None
    if args.per_base:
        pb_path = outdir / "coverage.perbase.bedgraph.gz"
        perbase_fh = gzip.open(pb_path, "wt")

    for chrom in chroms:
        base_depth, _, _ = calc_cov(cram, chrom, rf, per_base=True)
        length = cram.get_reference_length(chrom)
        hist = analysis.depth_histogram(base_depth)
        per_chrom_stats[chrom] = analysis.stats_from_hist(hist)
        lengths[chrom] = length

        if analysis.is_autosome(chrom):
            autosomal_hist = hist if autosomal_hist is None else autosomal_hist + hist

        # Stratified reduction: mask the per-base vector by each callability tier.
        easy_mask = None
        for label, beds in strata_beds.items():
            if chrom not in beds:
                continue
            mask = strata_mod.stratum_mask(length, *beds[chrom])
            h = analysis.depth_histogram(base_depth[mask])
            strata_hist[label] = h if label not in strata_hist else strata_hist[label] + h
            strata_bp[label] += int(mask.sum())
            if label == "easy":
                easy_mask = mask
                if analysis.is_autosome(chrom):
                    easy_autosomal_hist = h if easy_autosomal_hist is None else easy_autosomal_hist + h

        starts, ends, means = analysis.windowed_means(base_depth, args.window)
        # Per-window easy fraction (same reduceat trick on the mask) so the scatter
        # can drop repeat/centromere windows. No strata -> everything counts as easy.
        if easy_mask is not None:
            ef = np.add.reduceat(easy_mask.astype(np.int64), starts) / (ends - starts)
        else:
            ef = np.ones_like(means)
        for s, e, m, f in zip(starts.tolist(), ends.tolist(), means.tolist(), ef.tolist()):
            win_rows.append({"chrom": chrom, "start": s, "end": e, "mean": m, "easy_frac": f})

        if perbase_fh is not None:
            for s, e, d in analysis.rle_intervals(base_depth, skip_zero=True):
                perbase_fh.write(f"{chrom}\t{s}\t{e}\t{d}\n")

        del base_depth  # free the big vector before the next chromosome

    if perbase_fh is not None:
        perbase_fh.close()

    # Baseline = median depth of the diploid reference. Prefer the callability-
    # masked easy-autosomal positions (drops the segdup/centromere pileups that
    # otherwise inflate it); fall back to all autosomes, then to a length-weighted
    # mean if the subset has no autosomes.
    baseline = analysis.quantile_from_hist(easy_autosomal_hist, 0.5) if easy_autosomal_hist is not None else 0.0
    baseline_source = "easy-autosomal median"
    if not baseline and autosomal_hist is not None:
        baseline, baseline_source = analysis.quantile_from_hist(autosomal_hist, 0.5), "autosomal median"
    if not baseline:
        baseline = sum(per_chrom_stats[c].mean * lengths[c] for c in chroms) / sum(lengths.values()) or 1.0
        baseline_source = "length-weighted mean (no autosomes in subset)"

    # Attach CN + focal flag to windows and assemble the stats table.
    for w in win_rows:
        w["cn"] = analysis.copy_number(w["mean"], baseline)
        w["flag"] = qc.window_flag(w["cn"])

    stats_rows = []
    flagged: list[tuple[str, list[str]]] = []
    for chrom in chroms:
        s = per_chrom_stats[chrom]
        cn = analysis.copy_number(s.mean, baseline)
        flags = qc.chrom_flags(chrom, s, cn, baseline)
        if flags:
            flagged.append((chrom, flags))
        stats_rows.append({
            "chrom": chrom, "length": lengths[chrom],
            "mean": round(s.mean, 2), "median": round(s.median, 2),
            "sd": round(s.sd, 2), "cv": round(s.cv, 3),
            "mad": round(s.mad, 2), "robust_cv": round(s.robust_cv, 3),
            "q25": round(s.q25, 2), "q75": round(s.q75, 2), "iqr": round(s.iqr, 2),
            "breadth_1x": round(s.breadth[1], 4), "breadth_10x": round(s.breadth[10], 4),
            "breadth_20x": round(s.breadth[20], 4),
            "copy_number": round(cn, 2),
            "flags": ";".join(flags) if flags else "OK",
        })

    _write_stats(stats_rows, outdir / "coverage.stats.tsv")
    _write_windows(win_rows, outdir / "coverage.windows.bed")
    if strata_hist:
        _write_strata(strata_hist, strata_bp, outdir / "coverage.strata.tsv")

    chrom_means = {c: per_chrom_stats[c].mean for c in chroms}
    plots.bar_by_chromosome(chrom_means, outdir / "coverage.bar.png", baseline=baseline)
    # Mask the scatter to callable windows when strata are supplied (drops the
    # centromere/segdup spikes); otherwise plot every window.
    min_easy = 0.5 if "easy" in strata_hist else 0.0
    plots.scatter_windows(win_rows, outdir / "coverage.scatter.png", baseline=baseline,
                          min_easy_frac=min_easy)

    print(f"diploid baseline ({baseline_source}): {baseline:.2f}x")
    if flagged:
        print("abnormalities flagged:")
        for chrom, fl in flagged:
            print(f"  {chrom}: {';'.join(fl)}")
    else:
        print("no chromosome-level abnormalities flagged")
    n_focal = sum(1 for w in win_rows if w.get("flag", ".") != ".")
    print(f"focal windows flagged: {n_focal} / {len(win_rows)}")
    print(f"wrote stats/windows/plots to {outdir}/")
    if strata_hist:
        print(f"wrote callability strata table to {outdir}/coverage.strata.tsv")
    if args.per_base:
        print(f"wrote per-base RLE bedgraph to {outdir}/coverage.perbase.bedgraph.gz")


if __name__ == "__main__":
    main()
