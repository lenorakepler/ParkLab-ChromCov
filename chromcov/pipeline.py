"""
The single pipeline orchestrator: source -> reduce -> accumulate.

One entry point (`run`) replaces the two old orchestrators (`coverage.run_coverage`
and `qc_report.QCReport`). A run is parameterized by a `Depth` (MEAN vs FULL) and
a `Source` (the CRAM alignment, or existing on-disk tracks):

  Depth.MEAN  + Source.ALIGNMENT : per-contig width totals only, in memory, NO track
                                 written (`per_base=False`). The mean path.
  Depth.FULL  + Source.ALIGNMENT : per-contig per-base vector -> RLE track on disk
                                 (resumable, no-op if present), then a SEPARATE
                                 pass re-reads the tracks and reduces. The --full path.
  Depth.FULL  + Source.TRACKS  : reduce existing tracks with NO CRAM (lengths from
                                 the reference .fai). The plot/gather path.

Why the two alignment paths do not share one worker: a `pysam.AlignmentFile` cannot
be pickled across a `ProcessPoolExecutor` boundary, so every parallel worker opens
its OWN handle. And the FULL compute->reduce is decoupled THROUGH the on-disk track
on purpose -- that is what makes --full resumable and lets gather run CRAM-free. The
MEAN path stays track-free: a mean run must never write bedgraphs it wasn't asked for.
"""
from __future__ import annotations

import enum
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from . import policy, preflight
from .filtering import ReadFilter
from .io import alignment, track
from .kernel import calc_cov
from .categories import Strata
from .reduce import ChromDepth
from .result import RunResult


class Depth(enum.Enum):
    MEAN = "mean"   # widths only, no vector, no track
    FULL = "full"   # vector -> track -> reduce

class Source(enum.Enum):
    ALIGNMENT = "alignment"   # read from the CRAM
    TRACKS = "tracks"         # reduce existing on-disk tracks (no CRAM)

# --- worker units (module-level so they pickle for ProcessPoolExecutor) -------
# Each worker OPENS ITS OWN CRAM handle: a pysam.AlignmentFile is not picklable.

def _mean_partition(args):
    """
    MEAN worker: one contig's aligned-base total (per_base=False, no vector).
    """
    cfg, chrom = args
    with alignment.open_alignment(cfg) as reader:
        _, total_depth, _ = calc_cov(reader, chrom, ReadFilter.from_config(cfg), per_base=False)
        length = reader.get_reference_length(chrom)
    return chrom, length, int(total_depth)

def compute_partition(cfg, chrom, bedgraph_dir, force: bool = False) -> Path:
    """
    FULL compute unit: one contig's per-base depth -> RLE track. No-op (resume)
    if the track already exists and not force. Opens its own AlignmentFile so it is
    safe in a worker process.
    """
    if not force and track.has_bedgraph(bedgraph_dir, chrom):
        return track.bedgraph_path(bedgraph_dir, chrom)
    with alignment.open_alignment(cfg) as reader:
        base_depth, _, _ = calc_cov(reader, chrom, ReadFilter.from_config(cfg), per_base=True)
    return track.write_bedgraph(bedgraph_dir, chrom, ChromDepth(base_depth))

def _compute_worker(args) -> str:
    cfg, chrom, bedgraph_dir, force = args
    compute_partition(cfg, chrom, bedgraph_dir, force=force)
    return chrom

# --- the pipeline -------------------------------------------------------------

def run(cfg, *, depth: Depth = Depth.MEAN, source: Source = Source.ALIGNMENT,
        jobs: int = 1, force: bool = False, bedgraph_dir=None,
        categories: Strata | None = None, chroms: list[str] | None = None,
        skip_preflight: bool = False) -> RunResult:
    """
    Run one coverage pass and return the accumulated RunResult.

    Kept free of argparse/IO-formatting so tests/workflows can call it directly.
    `jobs > 1` computes contigs in parallel (each worker opens its own handle).
    """
    cats = categories if categories is not None else Strata.from_arg(cfg.strata)
    result = RunResult(cfg=cfg, depth=depth, categories=cats)

    # Resolve the contigs to process + their lengths, from the chosen source.
    if source is Source.ALIGNMENT:
        if not skip_preflight:
            result.preflight = preflight.preflight(cfg)
        with alignment.open_alignment(cfg) as reader:
            contig_list = alignment.list_contigs(reader, cfg)
            lengths = alignment.contig_lengths(reader, contig_list)

    else:  # Source.TRACKS -- CRAM-free
        contig_list = list(chroms) if chroms else track.bedgraph_chroms(bedgraph_dir)
        lengths = alignment.lengths_from_reference(cfg, contig_list)

    if depth is Depth.MEAN:
        _run_mean(cfg, result, contig_list, lengths, jobs)
        return result

    # FULL: (ALIGNMENT) ensure tracks exist, then reduce tracks in a separate pass.
    if source is Source.ALIGNMENT:
        _ensure_tracks(cfg, contig_list, bedgraph_dir, jobs, force)
    
    _reduce_tracks(result, contig_list, lengths, bedgraph_dir)
    policy.finalize(result)
    return result

def _run_mean(cfg, result, contigs, lengths, jobs) -> None:
    """
    MEAN source: per-contig aligned-base totals, no tracks, no vector.
    """
    if jobs and jobs > 1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for chrom, length, bases in ex.map(_mean_partition, [(cfg, c) for c in contigs]):
                result.add_mean(chrom, length, bases)
        return
    
    # Serial: reuse a single open handle across contigs (cheaper for many contigs).
    rf = ReadFilter.from_config(cfg)
    with alignment.open_alignment(cfg) as reader:
        for chrom in contigs:
            _, total_depth, _ = calc_cov(reader, chrom, rf, per_base=False)
            result.add_mean(chrom, lengths[chrom], int(total_depth))

def _ensure_tracks(cfg, contigs, bedgraph_dir, jobs, force) -> None:
    """
    Write per-base RLE tracks for every contig missing one (or all, if force).
    Resumable: an existing track is not recomputed.
    """
    Path(bedgraph_dir).mkdir(parents=True, exist_ok=True)
    todo = [c for c in contigs if force or not track.has_bedgraph(bedgraph_dir, c)]
    if not todo:
        return
    if jobs and jobs > 1:
        args = [(cfg, c, str(bedgraph_dir), force) for c in todo]
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            list(ex.map(_compute_worker, args))
    else:
        for c in todo:
            compute_partition(cfg, c, bedgraph_dir, force=force)

def _reduce_tracks(result, contigs, lengths, bedgraph_dir) -> None:
    """
    Separate reduce pass: re-read each track and fold it into the result. Serial
    in the main process (matches the in-place accumulation the reductions expect).
    """
    for chrom in (c for c in contigs if c in lengths):
        base_depth = track.read_bedgraph(track.bedgraph_path(bedgraph_dir, chrom), lengths[chrom])
        result.reduce_chrom(chrom, lengths[chrom], base_depth)
