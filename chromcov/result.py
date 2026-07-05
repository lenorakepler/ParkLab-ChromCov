"""
Normalized coverage result.

Both backends emit `list[ChromCoverage]`, so downstream code (TSV writer,
report, "diff the two backends" test) never has to know which tool produced the
numbers. `mean` uses the one definition both tools share: bases / length.

`base_depth` is the optional per-base int32 vector (only populated by the native
backend with per_base=True; mosdepth's per-base output is a BED we don't parse
here). Keeping it Optional means the summary path stays cheap.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np


@dataclass
class ChromCoverage:
    chrom: str
    length: int
    bases: int                              # total aligned bases (sum of per-base depth)
    backend: str
    base_depth: np.ndarray | None = field(default=None, repr=False)

    @property
    def mean(self) -> float:
        # Same definition in both tools; guard zero-length contigs.
        return self.bases / self.length if self.length else 0.0

    def as_row(self) -> dict:
        return {
            "chrom": self.chrom,
            "length": self.length,
            "bases": self.bases,
            "mean_coverage": round(self.mean, 2),   # mosdepth reports 2 decimals; match it
            "backend": self.backend,
        }


TSV_COLUMNS = ["chrom", "length", "bases", "mean_coverage", "backend"]


def write_tsv(rows: list[ChromCoverage], path) -> None:
    import csv
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r.as_row())
