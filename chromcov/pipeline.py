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

from . import policy, preflight, reduce_cache
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

# ==============================================================================
#  Dispatch to different pipelines to calculate or read previously calculated
#  coverage info from specified chromosomes and return the accumulated RunResult.
# ==============================================================================
def _noop_progress(done: int, total: int, chrom: str, phase: str,
                   nbases: int = 0, total_bases: int = 0) -> None:
    """Default progress hook: does nothing (silent, matches prior behaviour)."""

def run(cfg, *, depth: Depth = Depth.MEAN, source: Source = Source.ALIGNMENT,
        jobs: int = 1, force: bool = False, bedgraph_dir=None, cache_dir=None,
        categories: Strata | None = None, chroms: list[str] | None = None,
        skip_preflight: bool = False, progress=None) -> RunResult:
    """
    Calculate or fetch coverage for selected chromosomes and return a pooled RunResult.

    Source = ALIGNMENT -> Depth.MEAN -> Calculate ONLY mean per chrom, return.
                       -> Depth.FULL -> Calculate mean per base, save as bedgraph

    Source = TRACKS    ->               Get per base coverage info from bedgraphs


    Kept free of argparse/IO-formatting so tests/workflows can call it directly.
    `jobs > 1` computes contigs in parallel (each worker opens its own handle).
    `progress`, if given, is called `progress(done, total, chrom, phase)` after each
    contig -- the CLI uses it to render a status line; callers that want silence omit it.
    """
    report = progress or _noop_progress

    # Get stratification categories
    cats = categories if categories is not None else Strata.from_arg(cfg.strata)
    
    # Instantiate RunResult
    result = RunResult(cfg=cfg, depth=depth, categories=cats)

    # Get contig list and lengths from source
    # If source is CRAM alignment, do preflight validation
    if source is Source.ALIGNMENT:
        if not skip_preflight:
            result.preflight = preflight.preflight(cfg)
        with alignment.open_alignment(cfg) as reader:
            contig_list = alignment.list_contigs(reader, cfg)
            lengths = alignment.contig_lengths(reader, contig_list)

    # Otherwise, get contig list and lengths from bedgraph
    else:
        contig_list = list(chroms) if chroms else track.bedgraph_chroms(bedgraph_dir)
        lengths = alignment.lengths_from_reference(cfg, contig_list)

    # If just getting mean chrom. depth, calc and return.
    if depth is Depth.MEAN:
        _run_mean(cfg, result, contig_list, lengths, jobs, report)
        return result

    # If doing full, get + save per-base coverage for all chromosomes that do not
    # currently have an output bedgraph file
    if source is Source.ALIGNMENT:
        _ensure_tracks(cfg, contig_list, lengths, bedgraph_dir, jobs, force, report)

    #
    _reduce_tracks(result, contig_list, lengths, bedgraph_dir, report,
                   cache_dir=cache_dir, force=force)
    policy.finalize(result)
    return result

def _run_mean(cfg, result, contigs, lengths, jobs, report) -> None:
    """
    MEAN source: per-contig aligned-base totals, no tracks, no vector.
    """
    total = len(contigs)
    total_bp = sum(lengths.get(c, 0) for c in contigs)
    report(0, total, "", "scan", total_bases=total_bp)
    if jobs and jobs > 1:
        # Parallel: results arrive as workers finish, so report on completion.
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for done, (chrom, length, bases) in enumerate(
                    ex.map(_mean_partition, [(cfg, c) for c in contigs]), 1):
                result.add_mean(chrom, length, bases)
                report(done, total, chrom, "scan", nbases=length)
        return

    # Serial: reuse a single open handle across contigs (cheaper for many contigs).
    # Announce each contig BEFORE scanning it -- a single contig can take a while,
    # so reporting on completion would look like a stall.
    rf = ReadFilter.from_config(cfg)
    with alignment.open_alignment(cfg) as reader:
        for i, chrom in enumerate(contigs, 1):
            report(i, total, chrom, "scan", nbases=lengths.get(chrom, 0))
            _, total_depth, _ = calc_cov(reader, chrom, rf, per_base=False)
            result.add_mean(chrom, lengths[chrom], int(total_depth))

def _ensure_tracks(cfg, contigs, lengths, bedgraph_dir, jobs, force, report) -> None:
    """
    Write per-base RLE tracks for every contig missing one (or all, if force).
    Resumable: an existing track is not recomputed.

    Contigs are computed SMALLEST-FIRST: per-base depth for a whole contig is the
    slow step, so finishing a tiny contig early gives the caller a bp/sec rate (and
    thus a usable ETA) within seconds instead of after all of chr1. This reorders
    only the on-disk compute; the table is built later in `_reduce_tracks`, which
    iterates the natural contig order -- so output ordering is unaffected.
    """
    Path(bedgraph_dir).mkdir(parents=True, exist_ok=True)
    todo = [c for c in contigs if force or not track.has_bedgraph(bedgraph_dir, c)]
    todo.sort(key=lambda c: lengths.get(c, 0))
    total = len(todo)
    if not total:
        return
    total_bp = sum(lengths.get(c, 0) for c in todo)
    report(0, total, "", "depth", total_bases=total_bp)
    if jobs and jobs > 1:
        # Parallel: report as each worker's track lands (submission order == result
        # order, so smallest-first still yields the first ETA fastest).
        args = [(cfg, c, str(bedgraph_dir), force) for c in todo]
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            for done, chrom in enumerate(ex.map(_compute_worker, args), 1):
                report(done, total, chrom, "depth", nbases=lengths.get(chrom, 0))
    else:
        # Serial: announce each contig BEFORE computing it -- otherwise a large
        # contig looks like a hang.
        for i, c in enumerate(todo, 1):
            report(i, total, c, "depth", nbases=lengths.get(c, 0))
            compute_partition(cfg, c, bedgraph_dir, force=force)

def _reduce_tracks(result, contigs, lengths, bedgraph_dir, report,
                   cache_dir=None, force=False) -> None:
    """
    Separate reduce pass: fold each contig into the result. Serial in the main
    process (matches the in-place accumulation the reductions expect).

    Cache-aware: a contig whose reduced intermediate is cached (and matches the
    current reduce config) is loaded instead of re-reading and re-reducing its
    per-base track -- so re-plotting an existing run does no reduce work, and a run
    that adds contigs only reduces the new ones. Freshly reduced contigs are cached.
    """
    items = [c for c in contigs if c in lengths]
    total = len(items)
    total_bp = sum(lengths[c] for c in items)
    report(0, total, "", "reduce", total_bases=total_bp)
    key = reduce_cache.cache_key(result.cfg, result.categories) if cache_dir else None
    for i, chrom in enumerate(items, 1):
        report(i, total, chrom, "reduce", nbases=lengths[chrom])
        rc = None
        if cache_dir and not force:
            rc = reduce_cache.load(cache_dir, chrom, key)
        if rc is None:
            base_depth = track.read_bedgraph(track.bedgraph_path(bedgraph_dir, chrom), lengths[chrom])
            rc = result._reduce_one(chrom, lengths[chrom], base_depth)
            if cache_dir:
                reduce_cache.save(cache_dir, rc, key)
        result.fold(rc)
