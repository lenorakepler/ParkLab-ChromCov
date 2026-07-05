"""
Region strata for callability-aware coverage, keyed to Park Lab's
SMaHT_Regional_Categorization (GRCh38):

  easy       1000G strict accessibility mask (most reliably callable)
  difficult  PanMask pm151 minus 1000G mask (moderately mappable, artifact-prone)
  extreme    outside both masks (repetitive / structurally complex)

    raw.githubusercontent.com/parklab/SMaHT_Regional_Categorization/main/SMaHT_easy_hg38.bed.gz
    (+ SMaHT_difficult_hg38.bed.gz, SMaHT_extreme_hg38.bed.gz)

`Strata` holds one BED per label, loaded to per-chromosome interval arrays, and
turns those into a boolean position mask with the SAME finite-difference trick
the coverage calc uses (+1 at starts, -1 at ends, cumsum > 0) -- vectorized, and
it handles overlapping/unsorted intervals for free.
"""
from __future__ import annotations

import gzip
from pathlib import Path
import numpy as np

# Conventional label order, best -> worst callability.
STRATUM_ORDER = ("easy", "difficult", "extreme")


def _load_bed(path: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """BED[.gz] -> {chrom: (starts, ends)} as int arrays (0-based half-open)."""
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    starts: dict[str, list[int]] = {}
    ends: dict[str, list[int]] = {}
    with opener(path, "rt") as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            starts.setdefault(c, []).append(int(s))
            ends.setdefault(c, []).append(int(e))
    return {c: (np.asarray(starts[c]), np.asarray(ends[c])) for c in starts}


class Strata:
    """Callability tiers: label -> {chrom: (starts, ends)}, with position masks."""

    def __init__(self, beds: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]]):
        self._beds = beds

    @classmethod
    def from_arg(cls, spec: str | dict[str, str]) -> "Strata":
        """Build from 'easy=a.bed.gz,difficult=b.bed.gz' (order preserved) or a
        {label: path} mapping. Loads each BED to interval arrays."""
        if isinstance(spec, str):
            mapping: dict[str, str] = {}
            for item in spec.split(","):
                if not item:
                    continue
                label, _, path = item.partition("=")
                mapping[label.strip()] = path.strip()
        else:
            mapping = dict(spec)
        return cls({label: _load_bed(path) for label, path in mapping.items()})

    def labels(self) -> list[str]:
        return list(self._beds)

    def __bool__(self) -> bool:
        return bool(self._beds)

    def __contains__(self, label: str) -> bool:
        return label in self._beds

    def has(self, label: str, chrom: str) -> bool:
        return label in self._beds and chrom in self._beds[label]

    def mask(self, label: str, chrom: str, length: int) -> np.ndarray | None:
        """Boolean length-array: True where a position falls in any interval of
        `label` on `chrom`. None if this stratum has no intervals for the chrom."""
        if not self.has(label, chrom):
            return None
        starts, ends = self._beds[label][chrom]
        delta = np.zeros(length + 1, dtype=np.int32)
        np.add.at(delta, np.clip(starts, 0, length), 1)
        np.add.at(delta, np.clip(ends, 0, length), -1)
        return np.cumsum(delta[:-1]) > 0
