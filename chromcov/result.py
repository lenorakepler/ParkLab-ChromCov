"""
Normalized per-chromosome coverage result.

`ChromCoverage` is the combined record: the headline `mean` plus, when a per-base
pass was done, the full depth `stats` (median/sd/MAD/breadth/...) and the
analysis-level `copy_number` + `flags`. So there is ONE coverage table, not a
separate coverage.tsv and stats.tsv -- `--fast` just leaves the extra fields
unset (mean only). `write_tsv` picks its columns from whatever the rows carry.

`base_depth` is the optional per-base int32 vector (kept only transiently);
`mean = bases / length` (the definition the mosdepth add-on also uses).
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field

import numpy as np

from .analysis import ChromStats

BASE_COLUMNS = ["chrom", "length", "bases", "mean_coverage"]
STAT_COLUMNS = ["median", "sd", "cv", "mad", "robust_cv", "q25", "q75", "iqr",
                "breadth_1x", "breadth_10x", "breadth_20x"]
CN_COLUMNS = ["copy_number", "flags"]


@dataclass
class ChromCoverage:
    chrom: str
    length: int
    bases: int                                   # total aligned bases (sum of per-base depth)
    base_depth: np.ndarray | None = field(default=None, repr=False)
    stats: ChromStats | None = None              # per-base depth stats (unset in --fast)
    copy_number: float | None = None             # analysis-level (needs a baseline)
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


def columns(rows: list[ChromCoverage]) -> list[str]:
    """The columns present across `rows`: base always; stats + CN when carried."""
    cols = list(BASE_COLUMNS)
    if any(r.stats is not None for r in rows):
        cols += STAT_COLUMNS
    if any(r.copy_number is not None for r in rows):
        cols += CN_COLUMNS
    return cols


def write_rows(rows: list[ChromCoverage], fh) -> None:
    cols = columns(rows)
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    for r in rows:
        cells = r.as_row()
        w.writerow({c: cells.get(c, "") for c in cols})


def write_tsv(rows: list[ChromCoverage], path) -> None:
    with open(path, "w", newline="") as fh:
        write_rows(rows, fh)
