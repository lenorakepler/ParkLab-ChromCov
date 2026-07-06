"""

"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pysam

from . import perbase, plots, report, validate
from .calc_cov import calc_cov
from .config import Config, ReadFilter
from .depth import ChromDepth, DepthHistogram
from .qc_flags import chrom_flags, copy_number, is_autosome
from .report import ChromCoverage
from .strata import STRATUM_ORDER, Strata


def _read_filter(cfg: Config) -> ReadFilter:
    return ReadFilter(
        include_flags=cfg.include_flags,
        exclude_flags=cfg.exclude_flags,
        exclude_all_flags=cfg.exclude_all_flags,
        min_mapping_quality=cfg.min_mapping_quality,
    )

def compute_chrom(cfg: Config, chrom: str, bedgraph_dir: str | Path,
                  force: bool = False) -> Path:
    """MAP unit: compute one chromosome's per-base depth from the CRAM and write
    its bedgraph. No-op (resume) if the bedgraph already exists and not force.
    Opens its own AlignmentFile so it is safe to run in a worker process."""
    path = perbase.bedgraph_path(bedgraph_dir, chrom)
    if not force and perbase.has_bedgraph(bedgraph_dir, chrom):
        return path
    cram = pysam.AlignmentFile(
        str(cfg.cram), "rc",
        reference_filename=str(cfg.reference), index_filename=str(cfg.index))
    try:
        base_depth, _, _ = calc_cov(cram, chrom, _read_filter(cfg), per_base=True)
    finally:
        cram.close()
    return perbase.write_bedgraph(bedgraph_dir, chrom, ChromDepth(base_depth))

def _compute_worker(args) -> str:
    cfg, chrom, bedgraph_dir, force = args
    compute_chrom(cfg, chrom, bedgraph_dir, force=force)
    return chrom

class QCReport:
    """One `--full` run: compute per-chrom bedgraphs, then reduce them."""

    def __init__(self, cfg: Config, strata: Strata | None = None):
        self.cfg = cfg
        self.strata = strata if strata is not None else Strata.from_arg(cfg.strata)

        self.chroms: list[str] = []
        self.per_chrom_stats: dict = {}
        self.lengths: dict[str, int] = {}
        self.bases: dict[str, int] = {}
        self.win_rows: list[dict] = []
        self.autosomal_hist: DepthHistogram | None = None
        self.strata_hist: dict[str, DepthHistogram] = {}
        self.strata_bp: dict[str, int] = {label: 0 for label in self.strata.labels()}
        self.easy_autosomal_hist: DepthHistogram | None = None

        self._rows: list[ChromCoverage] = []
        self._windows_df = None
        self.flagged: list[tuple[str, list[str]]] = []
        self._baseline: float | None = None
        self._baseline_source: str = ""

    # --- driving the pass --------------------------------------------------

    def run(self, bedgraph_dir: str | Path, jobs: int = 1, force: bool = False) -> dict:
        """Preflight, compute all contigs (parallel, resumable), then gather.
        Returns the preflight report."""
        report = validate.preflight(self.cfg)
        cram = pysam.AlignmentFile(
            str(self.cfg.cram), "rc",
            reference_filename=str(self.cfg.reference), index_filename=str(self.cfg.index))
        try:
            chrom_list = self.cfg.select_contigs(cram.references)
            lengths = {c: cram.get_reference_length(c) for c in chrom_list}
        finally:
            cram.close()

        self._compute_all(chrom_list, bedgraph_dir, jobs, force)
        self._reduce(chrom_list, lengths, bedgraph_dir)
        return report

    def gather(self, bedgraph_dir: str | Path, chroms: list[str] | None = None) -> None:
        """Reduce whatever bedgraphs are present, no CRAM required (the `plot`
        entry point). Lengths come from the reference .fai."""
        chrom_list = list(chroms) if chroms else perbase.bedgraph_chroms(bedgraph_dir)
        lengths = self._lengths_from_reference(chrom_list)
        self._reduce(chrom_list, lengths, bedgraph_dir)

    def _lengths_from_reference(self, chroms: list[str]) -> dict[str, int]:
        fa = pysam.FastaFile(str(self.cfg.reference))
        try:
            names = set(fa.references)
            return {c: fa.get_reference_length(c) for c in chroms if c in names}
        finally:
            fa.close()

    def _compute_all(self, chrom_list, bedgraph_dir, jobs, force) -> None:
        Path(bedgraph_dir).mkdir(parents=True, exist_ok=True)
        todo = [c for c in chrom_list
                if force or not perbase.has_bedgraph(bedgraph_dir, c)]
        if not todo:
            return
        if jobs and jobs > 1:
            args = [(self.cfg, c, str(bedgraph_dir), force) for c in todo]
            with ProcessPoolExecutor(max_workers=jobs) as ex:
                list(ex.map(_compute_worker, args))
        else:
            for c in todo:
                compute_chrom(self.cfg, c, bedgraph_dir, force=force)

    def _reduce(self, chrom_list, lengths, bedgraph_dir) -> None:
        self.chroms = [c for c in chrom_list if c in lengths]
        for chrom in self.chroms:
            self._reduce_chrom(chrom, lengths[chrom], bedgraph_dir)
        self.finalize()

    def _reduce_chrom(self, chrom: str, length: int, bedgraph_dir) -> None:
        """Read one bedgraph back, reduce to hist + windows (+ strata), accumulate
        the pooled histograms, then drop the per-base vector."""
        base_depth = perbase.read_bedgraph(perbase.bedgraph_path(bedgraph_dir, chrom), length)
        total_depth = int(base_depth.sum())
        depth = ChromDepth(base_depth, cap=self.cfg.hist_cap,
                           breadth_thresholds=self.cfg.breadth_thresholds)

        hist = depth.histogram()
        self.per_chrom_stats[chrom] = hist.stats()
        self.lengths[chrom] = length
        self.bases[chrom] = total_depth

        if is_autosome(chrom):
            self.autosomal_hist = hist if self.autosomal_hist is None else self.autosomal_hist + hist

        masks: dict[str, np.ndarray] = {}
        for label in self.strata.labels():
            mask = self.strata.mask(label, chrom, length)
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
        # Per-window fraction in each stratum, plus the dominant tier -- the scatter
        # colors every window by its tier (easy/difficult/extreme) instead of only
        # keeping callable ones. `easy_frac` is retained for the optional
        # callable-only view; with no strata a window is easy_frac=1, stratum="".
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

        del base_depth

    # --- reductions --------------------------------------------------------

    def baseline(self) -> tuple[float, str]:
        """Diploid (CN=2) reference depth. Prefer the callability-masked
        easy-autosomal median (drops segdup/centromere pileups); fall back to all
        autosomes, then to a length-weighted mean if the subset has no autosomes."""
        if self._baseline is not None:
            return self._baseline, self._baseline_source

        val, src = 0.0, ""
        if self.cfg.baseline == "easy-autosomal-median" and self.easy_autosomal_hist is not None:
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
        """Assemble the combined per-chromosome rows (mean + per-base stats + copy
        number + QC flags) and the windowed table (vectorized CN + focal flag)."""
        baseline, _ = self.baseline()
        ploidy = self.cfg.ploidy

        self._windows_df = report.windows_frame(self.win_rows, baseline, ploidy, self.cfg.qc)

        self._rows = []
        self.flagged = []
        for chrom in self.chroms:
            s = self.per_chrom_stats[chrom]
            cn = copy_number(s.mean, baseline, ploidy)
            fl = chrom_flags(chrom, s, cn, baseline, self.cfg.qc)
            if fl:
                self.flagged.append((chrom, fl))
            self._rows.append(ChromCoverage(
                chrom=chrom, length=self.lengths[chrom], bases=self.bases[chrom],
                stats=s, copy_number=cn, flags=fl))

    def coverage_rows(self) -> list[ChromCoverage]:
        return self._rows

    # --- writing outputs ---------------------------------------------------

    def write_outputs(self, outdir: str | Path) -> dict[str, Path]:
        """Write the combined coverage table, windows, (strata), and plots under
        `outdir`. The per-base bedgraphs are a separate output (see perbase)."""
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        written: dict[str, Path] = {}

        coverage_tsv = outdir / "coverage.tsv"
        report.write_table(report.coverage_frame(self.coverage_rows()), coverage_tsv)
        written["coverage"] = coverage_tsv

        windows_bed = outdir / "coverage.windows.bed"
        report.write_table(self._windows_df, windows_bed)
        written["windows"] = windows_bed

        if self.strata_hist:
            order = [s for s in STRATUM_ORDER if s in self.strata_hist] + \
                    [s for s in self.strata_hist if s not in STRATUM_ORDER]
            strata_tsv = outdir / "coverage.strata.tsv"
            report.write_table(report.strata_frame(self.strata_hist, self.strata_bp, order), strata_tsv)
            written["strata"] = strata_tsv

        if self.cfg.plots:
            baseline, _ = self.baseline()
            chrom_means = {c: self.per_chrom_stats[c].mean for c in self.chroms}
            bar = outdir / "coverage.bar.png"
            plots.bar_by_chromosome(chrom_means, bar, baseline=baseline)
            written["bar"] = bar
            min_easy = self.cfg.scatter_min_easy_frac if "easy" in self.strata_hist else 0.0
            scatter = outdir / "coverage.scatter.png"
            plots.scatter_windows(self.win_rows, scatter, baseline=baseline,
                                  ploidy=self.cfg.ploidy, min_easy_frac=min_easy)
            written["scatter"] = scatter

        return written

    def summary_lines(self) -> list[str]:
        """Human-readable run summary (baseline, flagged chroms, focal counts)."""
        baseline, source = self.baseline()
        lines = [f"diploid baseline ({source}): {baseline:.2f}x"]
        if self.flagged:
            lines.append("abnormalities flagged:")
            lines += [f"  {chrom}: {';'.join(fl)}" for chrom, fl in self.flagged]
        else:
            lines.append("no chromosome-level abnormalities flagged")
        n_focal = int((self._windows_df["flag"] != ".").sum()) if (
            self._windows_df is not None and self._windows_df.height) else 0
        lines.append(f"focal windows flagged: {n_focal} / {len(self.win_rows)}")
        return lines
