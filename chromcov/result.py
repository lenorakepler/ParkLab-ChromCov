"""
Normalized coverage result: one row per chromosome. `mean = bases / length` (the
definition the mosdepth add-on also uses, so their tables are comparable).

`base_depth` is the optional per-base int32 vector (populated with per_base=True);
keeping it Optional means the summary path stays cheap.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class ChromCoverage:
    chrom: str
    length: int
    bases: int                              # total aligned bases (sum of per-base depth)
    base_depth: np.ndarray | None = field(default=None, repr=False)

    @property
    def mean(self) -> float:
        return self.bases / self.length if self.length else 0.0

    def as_row(self) -> dict:
        return {
            "chrom": self.chrom,
            "length": self.length,
            "bases": self.bases,
            "mean_coverage": round(self.mean, 2),   # 2 dp, matching mosdepth's summary
        }


TSV_COLUMNS = ["chrom", "length", "bases", "mean_coverage"]


def write_tsv(rows: list[ChromCoverage], path) -> None:
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r.as_row())
