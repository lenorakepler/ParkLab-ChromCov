"""
Per-base depth tracks: the expensive per-base coverage, written once as a
first-class, interoperable OUTPUT (not a hidden cache) and reused thereafter.

Each chromosome's per-base depth is stored as a standard `.per-base.bedgraph.gz`
(chrom/start/end/depth, RLE over constant-depth runs) under

    <root>/<coverage-key>/chrN.per-base.bedgraph.gz
    <root>/<coverage-key>/coverage.json      # provenance + filter params + per-chrom summary

so any genome tool (bedtools, IGV/UCSC after bigWig conversion, ...) can consume
them, and a later `chromcov coverage` can re-derive stats/windows/strata/plots
*without touching the CRAM again*. The directory is keyed by input identity +
read-filter params, so different filters never collide, and it is populated
per-chromosome (incremental: a later run adds only the chromosomes it needs).

Reuse must beat recompute, so store/load are vectorized: writing uses the
finite-difference change-points (`np.diff`) + `np.savetxt`; loading rebuilds the
vector with the same +d/-d/cumsum trick the coverage calc uses.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from . import provenance
from .analysis import ChromDepth
from .config import AnalysisConfig, CoverageConfig
from .output import _stable_input_id, run_params

TRACK_SUFFIX = ".per-base.bedgraph.gz"
SIDECAR = "coverage.json"


class PerBaseStore:
    """A directory of per-chromosome per-base depth tracks + a provenance sidecar,
    keyed by (input identity + read-filter params)."""

    def __init__(self, root: str | Path, config: CoverageConfig):
        self.root = Path(root)
        self.config = config
        self.key = self._key()
        self.dir = self.root / self.key

    def _key(self, length: int = 12) -> str:
        payload = {"params": run_params(self.config), "inputs": _stable_input_id(self.config)}
        blob = json.dumps(payload, sort_keys=True).encode()
        return hashlib.sha256(blob).hexdigest()[:length]

    # --- track paths / presence -------------------------------------------

    def track_path(self, chrom: str) -> Path:
        return self.dir / f"{chrom}{TRACK_SUFFIX}"

    def has(self, chrom: str) -> bool:
        return self.track_path(chrom).exists()

    def exists(self) -> bool:
        return (self.dir / SIDECAR).exists()

    def chroms(self) -> list[str]:
        return sorted(p.name[: -len(TRACK_SUFFIX)] for p in self.dir.glob(f"*{TRACK_SUFFIX}"))

    # --- store / load ------------------------------------------------------

    def store(self, chrom: str, depth: ChromDepth) -> Path:
        """Write chrom's per-base depth as an RLE bedgraph.gz (atomic). Vectorized
        over the vector's change-points; zero-depth runs are omitted (implicit)."""
        self.dir.mkdir(parents=True, exist_ok=True)
        d = depth.base_depth
        if d is None:
            raise ValueError("PerBaseStore.store needs a per-base depth vector (per_base=True)")

        if d.size == 0:
            starts = ends = vals = np.empty(0, dtype=np.int64)
        else:
            bounds = np.concatenate(([0], np.flatnonzero(np.diff(d)) + 1, [d.size]))
            starts, ends = bounds[:-1], bounds[1:]
            vals = d[starts]
            nz = vals != 0
            starts, ends, vals = starts[nz], ends[nz], vals[nz]

        path = self.track_path(chrom)
        tmp = path.with_suffix(path.suffix + ".part")
        with gzip.open(tmp, "wt") as fh:
            if vals.size:
                np.savetxt(fh, np.column_stack([starts, ends, vals]),
                           fmt=f"{chrom}\t%d\t%d\t%d")
        tmp.replace(path)
        return path

    def load(self, chrom: str, length: int) -> np.ndarray:
        """Reconstruct the per-base int32 depth vector from the stored RLE track.
        +d at each run start, -d at its end, cumsum -- exact since the runs
        partition the chromosome (the omitted gaps are depth 0)."""
        with gzip.open(self.track_path(chrom), "rb") as fh:
            tok = fh.read().split()
        diff = np.zeros(length + 1, dtype=np.int64)
        if tok:
            starts = np.array(tok[1::4]).astype(np.int64)   # skip col 0 (chrom name)
            ends = np.array(tok[2::4]).astype(np.int64)
            depths = np.array(tok[3::4]).astype(np.int64)
            np.add.at(diff, starts, depths)
            np.add.at(diff, ends, -depths)
        return np.cumsum(diff[:-1]).astype(np.int32)

    # --- sidecar -----------------------------------------------------------

    def write_sidecar(self, summary: dict[str, dict]) -> Path:
        """Write coverage.json: full provenance + read-filter params + a per-chrom
        summary (length/bases/mean) so the coverage table is available from the
        tracks alone, no CRAM required."""
        record = provenance.build_provenance(
            params=run_params(self.config),
            cram=self.config.cram,
            crai=self.config.index,
            reference=self.config.reference,
            outputs=[self.track_path(c) for c in summary],
        )
        record["coverage_key"] = self.key
        record["config"] = self.config.model_dump(mode="json")   # the resolved coverage config
        record["chromosomes"] = summary
        path = self.dir / SIDECAR
        path.write_text(json.dumps(record, indent=2, sort_keys=True))
        return path

    def read_summary(self) -> dict[str, dict]:
        """Per-chrom {length, bases, mean} from the sidecar (empty if none)."""
        path = self.dir / SIDECAR
        if not path.exists():
            return {}
        return json.loads(path.read_text()).get("chromosomes", {})


# --- Level 2 keys: an analysis run derived from a per-base coverage-key ------

def analysis_key(coverage_key: str, acfg: AnalysisConfig, strata_labels, length: int = 10) -> str:
    """Deterministic id for an analysis run: the coverage-key it derives from +
    the analysis params (window, strata set, baseline, ...)."""
    payload = {
        "coverage_key": coverage_key,
        "window": acfg.window,
        "hist_cap": acfg.hist_cap,
        "breadth_thresholds": list(acfg.breadth_thresholds),
        "ploidy": acfg.ploidy,
        "baseline": acfg.baseline,
        "strata": sorted(strata_labels),
    }
    blob = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:length]


def analysis_slug(acfg: AnalysisConfig, strata_labels, key: str) -> str:
    """Human-scannable analysis dir name: window + strata tag + hash."""
    tag = "strata_" + "_".join(sorted(strata_labels)) if strata_labels else "nostrata"
    return f"w{acfg.window}-{tag}-{key}"


def build_track(config: CoverageConfig, chrom: str, store: PerBaseStore) -> dict:
    """Compute ONE chromosome's per-base depth from the CRAM and store it as a
    track. The Snakemake scatter unit. Deliberately writes only the track file
    (not the shared sidecar) so parallel per-chrom jobs don't race on it -- the
    gather step (`analyze --per-base`) stamps coverage.json. Returns the per-chrom
    summary and is a no-op if the track already exists."""
    if store.has(chrom):
        return store.read_summary().get(chrom, {})

    import pysam
    from .read_filter import ReadFilter, calc_cov

    rf = ReadFilter(
        include_flags=config.include_flags,
        exclude_flags=config.exclude_flags,
        exclude_all_flags=config.exclude_all_flags,
        min_mapping_quality=config.min_mapping_quality,
    )
    cram = pysam.AlignmentFile(
        str(config.cram), "rc",
        reference_filename=str(config.reference),
        index_filename=str(config.index),
    )
    try:
        length = cram.get_reference_length(chrom)
        base_depth, total_depth, _ = calc_cov(cram, chrom, rf, per_base=True)
    finally:
        cram.close()

    store.store(chrom, ChromDepth(base_depth))
    return {"length": int(length), "bases": int(total_depth),
            "mean": round(int(total_depth) / length, 4) if length else 0.0}
