"""
Assemble + write the coverage tables with polars (was report.py), plus the
top-level output orchestration for a --full run (`write_outputs`, `summary_lines`).

Tables here are pure derivations of a finalized RunResult; `ChromCoverage` is the
one normalized per-chromosome row.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from ..categories import STRATUM_ORDER
from ..config.schema import QCThresholds
from ..reduce import ChromStats, DepthHistogram
from . import plots

STRATA_COLUMNS = ["stratum", "region_bp", "pct_of_analyzed", "mean", "median", "sd", "cv",
                  "breadth_1x", "breadth_10x", "breadth_20x"]


@dataclass
class ChromCoverage:
    chrom: str
    length: int
    bases: int                                   # total aligned bases (sum of per-base depth)
    stats: ChromStats | None = None              # per-base depth stats (unset for mean-only)
    copy_number: float | None = None             # needs a baseline
    flags: list[str] | None = None

    @property
    def mean(self) -> float:
        return self.bases / self.length if self.length else 0.0

    def as_row(self) -> dict:
        row = {
            "chrom": self.chrom,
            "length": self.length,
            "bases": self.bases,
            "mean_coverage": round(self.mean, 2),
        }
        if self.stats is not None:
            s = self.stats
            row.update({
                "median": round(s.median, 2), "sd": round(s.sd, 2), "cv": round(s.cv, 3),
                "mad": round(s.mad, 2), "robust_cv": round(s.robust_cv, 3),
                "q25": round(s.q25, 2), "q75": round(s.q75, 2), "iqr": round(s.iqr, 2),
                "breadth_1x": round(s.breadth.get(1, 0.0), 4),
                "breadth_10x": round(s.breadth.get(10, 0.0), 4),
                "breadth_20x": round(s.breadth.get(20, 0.0), 4),
            })
        if self.copy_number is not None:
            row["copy_number"] = round(self.copy_number, 2)
            row["flags"] = ";".join(self.flags) if self.flags else "OK"
        return row


def coverage_frame(rows: list[ChromCoverage]) -> pl.DataFrame:
    """Per-chromosome table. Columns are exactly what the rows carry: base only
    for the mean-only path; + stats/copy-number/flags for a --full run."""
    return pl.DataFrame([r.as_row() for r in rows]) if rows else pl.DataFrame()


def windows_frame(win_rows: list[dict], baseline: float, ploidy: int,
                  thr: QCThresholds = QCThresholds()) -> pl.DataFrame:
    """Per-window table with vectorized copy number + focal flag (the derivation
    that used to be a per-row Python loop over ~300k windows)."""
    if not win_rows:
        return pl.DataFrame(schema={"chrom": pl.String, "start": pl.Int64, "end": pl.Int64,
                                    "mean_depth": pl.Float64, "copy_number": pl.Float64,
                                    "flag": pl.String})
    df = pl.DataFrame(win_rows)   # chrom, start, end, mean, easy_frac
    cn = (ploidy * pl.col("mean") / baseline) if baseline else pl.lit(0.0)
    df = df.with_columns(cn.alias("copy_number"))
    df = df.with_columns(
        pl.when(pl.col("copy_number") <= thr.depleted_cn).then(pl.lit("DEPLETED"))
        .when(pl.col("copy_number") >= thr.gain_cn).then(pl.lit("GAIN"))
        .when(pl.col("copy_number") <= thr.loss_cn).then(pl.lit("LOSS"))
        .otherwise(pl.lit(".")).alias("flag")
    )
    return df.select(
        pl.col("chrom"), pl.col("start"), pl.col("end"),
        pl.col("mean").round(2).alias("mean_depth"),
        pl.col("copy_number").round(3),
        pl.col("flag"),
    )


def strata_frame(strata_hist: dict[str, DepthHistogram], strata_bp: dict[str, int],
                 order: list[str]) -> pl.DataFrame:
    """Per-callability-tier coverage summary."""
    total_bp = sum(strata_bp.values()) or 1
    records = []
    for label in order:
        s = strata_hist[label].stats()
        records.append({
            "stratum": label, "region_bp": strata_bp[label],
            "pct_of_analyzed": round(100 * strata_bp[label] / total_bp, 2),
            "mean": round(s.mean, 2), "median": round(s.median, 2),
            "sd": round(s.sd, 2), "cv": round(s.cv, 3),
            "breadth_1x": round(s.breadth[1], 4), "breadth_10x": round(s.breadth[10], 4),
            "breadth_20x": round(s.breadth[20], 4),
        })
    return pl.DataFrame(records, schema=STRATA_COLUMNS) if records else pl.DataFrame()


def write_table(df: pl.DataFrame, dest: str | Path | None) -> None:
    """Write a frame as TSV. dest None or '-' -> stdout; else to the file path."""
    if dest in (None, "-"):
        sys.stdout.write(df.write_csv(separator="\t"))
        return
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    df.write_csv(str(dest), separator="\t")


def write_outputs(result, outdir: str | Path) -> dict[str, Path]:
    """Write the combined coverage table, windows, (strata), and plots for a
    finalized --full RunResult. The per-base tracks are a separate output."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    baseline, _ = result.baseline

    coverage_tsv = outdir / "coverage.tsv"
    write_table(coverage_frame(result.coverage_rows()), coverage_tsv)
    written["coverage"] = coverage_tsv

    windows_bed = outdir / "coverage.windows.bed"
    write_table(windows_frame(result.win_rows, baseline, result.cfg.ploidy, result.cfg.qc), windows_bed)
    written["windows"] = windows_bed

    if result.strata_hist:
        order = [s for s in STRATUM_ORDER if s in result.strata_hist] + \
                [s for s in result.strata_hist if s not in STRATUM_ORDER]
        strata_tsv = outdir / "coverage.strata.tsv"
        write_table(strata_frame(result.strata_hist, result.strata_bp, order), strata_tsv)
        written["strata"] = strata_tsv

    if result.cfg.plots:
        chrom_means = {c: result.per_chrom_stats[c].mean for c in result.chroms}
        bar = outdir / "coverage.bar.png"
        plots.bar_by_chromosome(chrom_means, bar, baseline=baseline)
        written["bar"] = bar
        min_easy = result.cfg.scatter_min_easy_frac if "easy" in result.strata_hist else 0.0
        scatter = outdir / "coverage.scatter.png"
        scatter_out = plots.scatter_windows(result.win_rows, scatter, baseline=baseline,
                                            ploidy=result.cfg.ploidy, min_easy_frac=min_easy,
                                            cap_cn=result.cfg.scatter_cap_cn)
        written["scatter"] = scatter_out["png"]
        written["scatter_html"] = scatter_out["html"]

    return written


def summary_lines(result) -> list[str]:
    """Human-readable run summary (baseline, flagged chroms, focal counts)."""
    baseline, source = result.baseline
    lines = [f"diploid baseline ({source}): {baseline:.2f}x"]
    if result.flagged:
        lines.append("abnormalities flagged:")
        lines += [f"  {chrom}: {';'.join(fl)}" for chrom, fl in result.flagged]
    else:
        lines.append("no chromosome-level abnormalities flagged")
    wdf = windows_frame(result.win_rows, baseline, result.cfg.ploidy, result.cfg.qc)
    n_focal = int((wdf["flag"] != ".").sum()) if wdf.height else 0
    lines.append(f"focal windows flagged: {n_focal} / {len(result.win_rows)}")
    return lines
