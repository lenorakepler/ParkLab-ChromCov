"""
Coverage analysis pipeline: per-chromosome stats + windowed track + copy number
+ callability strata + plots, all from one memory-bounded pass over the native
per-base depth.

`CoverageAnalysis` owns the accumulators that used to be loose locals in
`analyze.py::main()` (`per_chrom_stats`, pooled autosomal/strata histograms,
windowed rows, ...). Each chromosome is reduced to a few-KB histogram + a
windowed track via `ChromDepth`, then its per-base vector is dropped, so peak
memory stays ~one chromosome.

Because it already does the full per-base pass, `analyze` is a *superset* of
`coverage`: `write_outputs` also emits the plain per-chromosome coverage table
(`coverage.tsv`), so one whole-genome run yields the deliverable + the QC suite.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pysam

from . import plots, qc, validate
from .analysis import (
    ChromDepth,
    ChromStats,
    DepthHistogram,
    copy_number,
    is_autosome,
)
from .config import AnalysisConfig, CoverageConfig
from .read_filter import ReadFilter, calc_cov
from .result import ChromCoverage, write_tsv
from .strata import STRATUM_ORDER, Strata

STATS_COLUMNS = [
    "chrom", "length", "mean", "median", "sd", "cv", "mad", "robust_cv",
    "q25", "q75", "iqr", "breadth_1x", "breadth_10x", "breadth_20x",
    "copy_number", "flags",
]

STRATA_COLUMNS = ["stratum", "region_bp", "pct_of_analyzed", "mean", "median", "sd", "cv",
                  "breadth_1x", "breadth_10x", "breadth_20x"]


class CoverageAnalysis:
    """One memory-bounded pass over the CRAM -> stats/windows/strata/plots."""

    def __init__(self, config: CoverageConfig, analysis: AnalysisConfig | None = None,
                 strata: Strata | None = None):
        self.config = config
        self.acfg = analysis or AnalysisConfig()
        self.strata = strata if strata is not None else Strata.from_arg(self.acfg.strata)

        # --- accumulators (formerly loose locals in analyze.main) ---
        self.chroms: list[str] = []
        self.per_chrom_stats: dict[str, ChromStats] = {}
        self.lengths: dict[str, int] = {}
        self.bases: dict[str, int] = {}          # total aligned bases (for the coverage table)
        self.win_rows: list[dict] = []
        self.autosomal_hist: DepthHistogram | None = None
        self.strata_hist: dict[str, DepthHistogram] = {}
        self.strata_bp: dict[str, int] = {label: 0 for label in self.strata.labels()}
        self.easy_autosomal_hist: DepthHistogram | None = None

        # --- per-base track reuse bookkeeping ---
        self.track_hits: list[str] = []      # chroms read back from stored tracks
        self.track_misses: list[str] = []    # chroms computed fresh from the CRAM

        # --- results (filled by finalize) ---
        self.stats_rows: list[dict] = []
        self.flagged: list[tuple[str, list[str]]] = []
        self._baseline: float | None = None
        self._baseline_source: str = ""

    # --- driving the pass --------------------------------------------------

    def run(self, chroms: list[str] | None = None, store=None,
            write_tracks: bool = False) -> dict:
        """Preflight, loop chromosomes, finalize. Returns the preflight report.
        Whole genome (contig globs) unless `chroms` is given.

        Per-chromosome per-base depth comes from `store` (a PerBaseStore) when a
        track exists there, else it's computed from the CRAM; with `write_tracks`
        a freshly computed chromosome is written back as a track output. If every
        requested chromosome is already in the store, the CRAM is never opened
        (lengths come from the store sidecar), so analysis can run without it."""
        cram_available = Path(self.config.cram).exists()
        cram = rf = None
        lengths: dict[str, int] = {}

        if cram_available:
            report = validate.preflight(self.config)
            rf = ReadFilter(
                include_flags=self.config.include_flags,
                exclude_flags=self.config.exclude_flags,
                exclude_all_flags=self.config.exclude_all_flags,
                min_mapping_quality=self.config.min_mapping_quality,
            )
            cram = pysam.AlignmentFile(
                str(self.config.cram), "rc",
                reference_filename=str(self.config.reference),
                index_filename=str(self.config.index),
            )
            self.chroms = list(chroms) if chroms else self.config.select_contigs(cram.references)
            lengths = {c: cram.get_reference_length(c) for c in self.chroms}
        else:
            if store is None or not store.exists():
                raise RuntimeError(
                    "CRAM not found and no per-base tracks to read from: "
                    f"{self.config.cram}"
                )
            report = {"reference_check": {"status": "reused per-base tracks (no CRAM)"}}
            summary = store.read_summary()
            self.chroms = list(chroms) if chroms else store.chroms()
            lengths = {c: summary[c]["length"] for c in self.chroms if c in summary}

        try:
            for chrom in self.chroms:
                length = lengths.get(chrom) or (cram.get_reference_length(chrom) if cram else None)
                if length is None:
                    raise RuntimeError(f"cannot determine length of {chrom!r} (not in CRAM or tracks)")
                self.process_chrom(chrom, length, cram=cram, rf=rf,
                                   store=store, write_tracks=write_tracks)
        finally:
            if cram is not None:
                cram.close()

        # Refresh the track sidecar with whatever chromosomes we now have (merging
        # with any already recorded, so a subset run doesn't drop other chroms).
        if store is not None and write_tracks:
            merged = {**store.read_summary(), **{
                c: {"length": self.lengths[c], "bases": self.bases[c],
                    "mean": round(self.per_chrom_stats[c].mean, 4)}
                for c in self.chroms
            }}
            store.write_sidecar(merged)

        self.finalize()
        return report

    def process_chrom(self, chrom: str, length: int, cram=None, rf: ReadFilter | None = None,
                      store=None, write_tracks: bool = False) -> None:
        """Reduce one chromosome to hist + windows (+ strata), accumulate pooled
        histograms, then drop the per-base vector. Per-base depth is loaded from
        `store` when present, else computed from the CRAM."""
        from_store = store is not None and store.has(chrom)
        if from_store:
            base_depth = store.load(chrom, length)
            total_depth = int(base_depth.sum())
            self.track_hits.append(chrom)
        else:
            if cram is None:
                raise RuntimeError(
                    f"chromosome {chrom!r} is not in the per-base tracks and no CRAM "
                    "is available to compute it"
                )
            base_depth, total_depth, _ = calc_cov(cram, chrom, rf, per_base=True)
            total_depth = int(total_depth)
            self.track_misses.append(chrom)

        depth = ChromDepth(base_depth, cap=self.acfg.hist_cap,
                           breadth_thresholds=self.acfg.breadth_thresholds)

        hist = depth.histogram()
        self.per_chrom_stats[chrom] = hist.stats()
        self.lengths[chrom] = length
        self.bases[chrom] = int(total_depth)

        if is_autosome(chrom):
            self.autosomal_hist = hist if self.autosomal_hist is None else self.autosomal_hist + hist

        # Stratified reduction: mask the per-base vector by each callability tier.
        easy_mask = None
        for label in self.strata.labels():
            mask = self.strata.mask(label, chrom, length)
            if mask is None:
                continue
            h = depth.masked(mask).histogram()
            self.strata_hist[label] = h if label not in self.strata_hist else self.strata_hist[label] + h
            self.strata_bp[label] += int(mask.sum())
            if label == "easy":
                easy_mask = mask
                if is_autosome(chrom):
                    self.easy_autosomal_hist = h if self.easy_autosomal_hist is None else self.easy_autosomal_hist + h

        starts, ends, means = depth.windowed_means(self.acfg.window)
        # Per-window easy fraction (same reduceat trick on the mask) so the scatter
        # can drop repeat/centromere windows. No strata -> everything counts as easy.
        if easy_mask is not None:
            ef = np.add.reduceat(easy_mask.astype(np.int64), starts) / (ends - starts)
        else:
            ef = np.ones_like(means)
        for s, e, m, f in zip(starts.tolist(), ends.tolist(), means.tolist(), ef.tolist()):
            self.win_rows.append({"chrom": chrom, "start": s, "end": e, "mean": m, "easy_frac": f})

        # Persist the per-base depth as a track output (only when freshly computed
        # and asked to); a reused track is already on disk.
        if store is not None and write_tracks and not from_store:
            store.store(chrom, depth)

        del base_depth  # free the big vector before the next chromosome

    # --- reductions --------------------------------------------------------

    def baseline(self) -> tuple[float, str]:
        """Diploid (CN=2) reference depth. Prefer the callability-masked
        easy-autosomal median (drops segdup/centromere pileups); fall back to all
        autosomes, then to a length-weighted mean if the subset has no autosomes."""
        if self._baseline is not None:
            return self._baseline, self._baseline_source

        val, src = 0.0, ""
        if self.acfg.baseline == "easy-autosomal-median" and self.easy_autosomal_hist is not None:
            val, src = self.easy_autosomal_hist.quantile(0.5), "easy-autosomal median"
        if not val and self.autosomal_hist is not None:
            val, src = self.autosomal_hist.quantile(0.5), "autosomal median"
        if not val:
            total_len = sum(self.lengths.values()) or 1
            val = sum(self.per_chrom_stats[c].mean * self.lengths[c] for c in self.chroms) / total_len or 1.0
            src = "length-weighted mean (no autosomes in subset)"

        self._baseline, self._baseline_source = val, src
        return val, src

    def finalize(self) -> None:
        """Attach copy number + QC flags to windows and chromosomes."""
        baseline, _ = self.baseline()
        ploidy = self.acfg.ploidy

        for w in self.win_rows:
            w["cn"] = copy_number(w["mean"], baseline, ploidy)
            w["flag"] = qc.window_flag(w["cn"])

        self.stats_rows = []
        self.flagged = []
        for chrom in self.chroms:
            s = self.per_chrom_stats[chrom]
            cn = copy_number(s.mean, baseline, ploidy)
            flags = qc.chrom_flags(chrom, s, cn, baseline)
            if flags:
                self.flagged.append((chrom, flags))
            self.stats_rows.append({
                "chrom": chrom, "length": self.lengths[chrom],
                "mean": round(s.mean, 2), "median": round(s.median, 2),
                "sd": round(s.sd, 2), "cv": round(s.cv, 3),
                "mad": round(s.mad, 2), "robust_cv": round(s.robust_cv, 3),
                "q25": round(s.q25, 2), "q75": round(s.q75, 2), "iqr": round(s.iqr, 2),
                "breadth_1x": round(s.breadth[1], 4), "breadth_10x": round(s.breadth[10], 4),
                "breadth_20x": round(s.breadth[20], 4),
                "copy_number": round(cn, 2),
                "flags": ";".join(flags) if flags else "OK",
            })

    def coverage_rows(self) -> list[ChromCoverage]:
        """The plain per-chromosome coverage table (the `coverage` deliverable),
        derived from the same pass -- this is what makes `analyze` a superset."""
        return [
            ChromCoverage(chrom=c, length=self.lengths[c], bases=self.bases[c])
            for c in self.chroms
        ]

    # --- writing outputs ---------------------------------------------------

    def write_outputs(self, outdir: Path) -> dict[str, Path]:
        """Write the coverage table, stats, windows, (strata), and plots under
        `outdir` (a Level-2 analysis run). The per-base tracks are a separate
        Level-1 output (see PerBaseStore). Returns a map of name -> path."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        coverage_tsv = outdir / "coverage.tsv"
        write_tsv(self.coverage_rows(), coverage_tsv)
        written["coverage"] = coverage_tsv

        stats_tsv = outdir / "coverage.stats.tsv"
        self._write_stats(stats_tsv)
        written["stats"] = stats_tsv

        windows_bed = outdir / "coverage.windows.bed"
        self._write_windows(windows_bed)
        written["windows"] = windows_bed

        if self.strata_hist:
            strata_tsv = outdir / "coverage.strata.tsv"
            self._write_strata(strata_tsv)
            written["strata"] = strata_tsv

        if self.acfg.plots:
            baseline, _ = self.baseline()
            chrom_means = {c: self.per_chrom_stats[c].mean for c in self.chroms}
            bar = outdir / "coverage.bar.png"
            plots.bar_by_chromosome(chrom_means, bar, baseline=baseline)
            written["bar"] = bar
            # Mask the scatter to callable windows when strata are supplied.
            min_easy = self.acfg.scatter_min_easy_frac if "easy" in self.strata_hist else 0.0
            scatter = outdir / "coverage.scatter.png"
            plots.scatter_windows(self.win_rows, scatter, baseline=baseline,
                                  ploidy=self.acfg.ploidy, min_easy_frac=min_easy)
            written["scatter"] = scatter

        return written

    def _write_stats(self, path: Path) -> None:
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=STATS_COLUMNS, delimiter="\t")
            w.writeheader()
            w.writerows(self.stats_rows)

    def _write_windows(self, path: Path) -> None:
        with path.open("w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["chrom", "start", "end", "mean_depth", "copy_number", "flag"])
            for x in self.win_rows:
                w.writerow([x["chrom"], x["start"], x["end"], f"{x['mean']:.2f}",
                            f"{x['cn']:.3f}", x.get("flag", ".")])

    def _write_strata(self, path: Path) -> None:
        total_bp = sum(self.strata_bp.values()) or 1
        order = [s for s in STRATUM_ORDER if s in self.strata_hist] + \
                [s for s in self.strata_hist if s not in STRATUM_ORDER]
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=STRATA_COLUMNS, delimiter="\t")
            w.writeheader()
            for label in order:
                s = self.strata_hist[label].stats()
                w.writerow({
                    "stratum": label, "region_bp": self.strata_bp[label],
                    "pct_of_analyzed": round(100 * self.strata_bp[label] / total_bp, 2),
                    "mean": round(s.mean, 2), "median": round(s.median, 2),
                    "sd": round(s.sd, 2), "cv": round(s.cv, 3),
                    "breadth_1x": round(s.breadth[1], 4), "breadth_10x": round(s.breadth[10], 4),
                    "breadth_20x": round(s.breadth[20], 4),
                })

    def summary_lines(self) -> list[str]:
        """Human-readable run summary (baseline, flagged chroms, focal counts)."""
        baseline, source = self.baseline()
        lines = []
        if self.track_hits:
            lines.append(f"per-base tracks reused: {len(self.track_hits)}/{len(self.chroms)} "
                         f"chrom(s) ({len(self.track_misses)} computed from CRAM)")
        lines.append(f"diploid baseline ({source}): {baseline:.2f}x")
        if self.flagged:
            lines.append("abnormalities flagged:")
            lines += [f"  {chrom}: {';'.join(fl)}" for chrom, fl in self.flagged]
        else:
            lines.append("no chromosome-level abnormalities flagged")
        n_focal = sum(1 for w in self.win_rows if w.get("flag", ".") != ".")
        lines.append(f"focal windows flagged: {n_focal} / {len(self.win_rows)}")
        return lines
