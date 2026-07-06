"""
Record per-base coverage
"""
from __future__ import annotations

import gzip
import os
from pathlib import Path

import numpy as np

from .depth import ChromDepth

BEDGRAPH_SUFFIX = ".per-base.bedgraph.gz"

def bedgraph_path(directory: str | Path, chrom: str) -> Path:
    return Path(directory) / f"{chrom}{BEDGRAPH_SUFFIX}"

def has_bedgraph(directory: str | Path, chrom: str) -> bool:
    return bedgraph_path(directory, chrom).exists()

def bedgraph_chroms(directory: str | Path) -> list[str]:
    """Chromosomes with a bedgraph present -- what the gather/plot step reads."""
    directory = Path(directory)
    if not directory.exists():
        return []
    return sorted(p.name[: -len(BEDGRAPH_SUFFIX)]
                  for p in directory.glob(f"*{BEDGRAPH_SUFFIX}"))

def write_bedgraph(directory: str | Path, chrom: str, depth: ChromDepth) -> Path:
    """Write chrom's per-base depth as an RLE bedgraph.gz, atomically. Vectorized
    over the vector's change-points; zero-depth runs are omitted (implicit)."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    d = depth.base_depth
    if d is None:
        raise ValueError("write_bedgraph needs a per-base depth vector (per_base=True)")

    if d.size == 0:
        starts = ends = vals = np.empty(0, dtype=np.int64)
    else:
        bounds = np.concatenate(([0], np.flatnonzero(np.diff(d)) + 1, [d.size]))
        starts, ends = bounds[:-1], bounds[1:]
        vals = d[starts]
        nz = vals != 0
        starts, ends, vals = starts[nz], ends[nz], vals[nz]

    path = bedgraph_path(directory, chrom)
    tmp = path.with_suffix(path.suffix + ".part")
    with gzip.open(tmp, "wt") as fh:
        if vals.size:
            np.savetxt(fh, np.column_stack([starts, ends, vals]),
                       fmt=f"{chrom}\t%d\t%d\t%d")
    os.replace(tmp, path)   # atomic: a complete file appears in one step
    return path

def read_bedgraph(path: str | Path, length: int) -> np.ndarray:
    """Reconstruct the per-base int32 depth vector from a stored RLE bedgraph.

    +d at each run start, -d at its end, cumsum -- exact since the runs partition
    the chromosome (omitted gaps are depth 0). `length` is needed because the RLE
    drops trailing zeros; the caller supplies it from the reference (.fai)/CRAM.
    """
    with gzip.open(path, "rb") as fh:
        tok = fh.read().split()
    diff = np.zeros(length + 1, dtype=np.int64)
    if tok:
        starts = np.array(tok[1::4]).astype(np.int64)    # skip col 0 (chrom name)
        ends = np.array(tok[2::4]).astype(np.int64)
        depths = np.array(tok[3::4]).astype(np.int64)
        np.add.at(diff, starts, depths)
        np.add.at(diff, ends, -depths)
    return np.cumsum(diff[:-1]).astype(np.int32)
