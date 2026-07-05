"""
Region strata for callability-aware coverage, keyed to Park Lab's
SMaHT_Regional_Categorization (GRCh38):

  easy       1000G strict accessibility mask (most reliably callable)
  difficult  PanMask pm151 minus 1000G mask (moderately mappable, artifact-prone)
  extreme    outside both masks (repetitive / structurally complex)

    raw.githubusercontent.com/parklab/SMaHT_Regional_Categorization/main/SMaHT_easy_hg38.bed.gz
    (+ SMaHT_difficult_hg38.bed.gz, SMaHT_extreme_hg38.bed.gz)

Each stratum is one BED. We load it to per-chromosome interval arrays and turn
those into a boolean position mask with the SAME finite-difference trick the
coverage calc uses (+1 at starts, -1 at ends, cumsum > 0) -- vectorized, and it
handles overlapping/unsorted intervals for free.
"""
from __future__ import annotations

import gzip
from pathlib import Path
import numpy as np

# Conventional label order, best -> worst callability.
STRATUM_ORDER = ("easy", "difficult", "extreme")


def load_bed(path: str | Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
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


def parse_strata_arg(spec: str) -> dict[str, str]:
    """'easy=a.bed.gz,difficult=b.bed.gz' -> {label: path}, order preserved."""
    out: dict[str, str] = {}
    for item in spec.split(","):
        if not item:
            continue
        label, _, path = item.partition("=")
        out[label.strip()] = path.strip()
    return out


def stratum_mask(length: int, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    """Boolean length-array: True where a position falls in any interval."""
    delta = np.zeros(length + 1, dtype=np.int32)
    np.add.at(delta, np.clip(starts, 0, length), 1)
    np.add.at(delta, np.clip(ends, 0, length), -1)
    return np.cumsum(delta[:-1]) > 0
