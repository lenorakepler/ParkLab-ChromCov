"""
Level 3: compare analysis runs.

Walks the analysis runs under an output root (out/<coverage-key>/<analysis-key>/),
reads each run's combined `coverage.tsv` + `run.json`, and pivots a chosen metric
to a chrom x run table. Because analysis runs off the same coverage-key share the
per-base stats, the interesting comparison is usually `copy_number` (its baseline
shifts with `--strata`) -- so "stratified vs unstratified" sits side by side.

This is the analysis-run analogue of `RunStore.collate` (which compares the
`runs/<id>/` coverage-table archives from `coverage --write`).
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

from .output import karyotypic_key

# Numeric columns of the combined coverage.tsv that make sense to pivot.
METRICS = ("mean_coverage", "median", "sd", "cv", "mad", "robust_cv",
           "breadth_1x", "breadth_10x", "breadth_20x", "copy_number")


class AnalysisRun:
    """One analysis-run directory: its combined table + run.json params."""

    def __init__(self, run_dir: Path):
        self.dir = Path(run_dir)
        # run_id = <coverage-key>/<analysis-slug>, unique and informative.
        self.run_id = f"{self.dir.parent.name}/{self.dir.name}"
        sidecar = self.dir / "run.json"
        self.sidecar = json.loads(sidecar.read_text()) if sidecar.exists() else {}
        self.rows = self._read_table()

    def _read_table(self) -> dict[str, dict]:
        rows: dict[str, dict] = {}
        tsv = self.dir / "coverage.tsv"
        if not tsv.exists():
            return rows
        with tsv.open() as fh:
            reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
            for rec in reader:
                rows[rec["chrom"]] = rec
        return rows

    @property
    def params(self) -> dict:
        cfg = self.sidecar.get("config", {})
        cov = cfg.get("coverage", {})
        ana = cfg.get("analysis", {})
        # Prefer the *effective* baseline source (from the run's baseline record)
        # over the configured mode -- an unstratified run falls back from
        # easy-autosomal to plain autosomal median.
        baseline = (self.sidecar.get("baseline") or {}).get("source") or ana.get("baseline")
        return {
            "coverage_key": self.sidecar.get("coverage_key"),
            "min_mapq": cov.get("min_mapping_quality"),
            "window": ana.get("window"),
            "baseline": baseline,
            "strata": sorted((ana.get("strata") or {})),
        }


def find_runs(root: str | Path = "out") -> list[AnalysisRun]:
    """Every analysis run under `root` (out/<coverage-key>/<analysis-key>/run.json)."""
    root = Path(root)
    return [AnalysisRun(p.parent) for p in sorted(root.glob("*/*/run.json"))]


def pivot(runs: list[AnalysisRun], metric: str = "copy_number"):
    """Wide table: chrom -> {run_id: metric value}, karyotypically ordered.

    Returns (run_ids, table). Runs missing the metric (e.g. --fast archives that
    only have mean_coverage) simply contribute no cell.
    """
    run_ids = [r.run_id for r in runs]
    table: dict[str, dict[str, str]] = {}
    for r in runs:
        for chrom, rec in r.rows.items():
            val = rec.get(metric, "")
            if val != "":
                table.setdefault(chrom, {})[r.run_id] = val
    ordered = {c: table[c] for c in sorted(table, key=karyotypic_key)}
    return run_ids, ordered
