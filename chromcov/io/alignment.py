"""
Alignment (CRAM) access: a context-managed reader, plus contig listing
and length lookup.

NOTE on "open once": a serial pass opens ONE handle and reuses it (via
`open_alignment`). Parallel workers CANNOT share a handle -- a `pysam.AlignmentFile`
wraps an open htslib handle and is not picklable across a process boundary -- so
each worker opens its own (see pipeline). This module provides the handle; who
holds it (main process vs worker) is the caller's decision.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator

import pysam


def open_reader(cfg) -> pysam.AlignmentFile:
    """Open a CRAM reader for this run's inputs. The caller owns closing it (use
    `open_alignment` for a context-managed handle)."""
    return pysam.AlignmentFile(
        str(cfg.cram), "rc",
        reference_filename=str(cfg.reference), index_filename=str(cfg.index))


@contextmanager
def open_alignment(cfg) -> Iterator[pysam.AlignmentFile]:
    """Context-managed CRAM reader, closed on exit. Safe to reuse across contigs
    within a single (serial) pass."""
    reader = open_reader(cfg)
    try:
        yield reader
    finally:
        reader.close()


def read_header(cfg) -> dict:
    """The CRAM header as a dict (for preflight: sort order, @SQ M5 tags)."""
    with open_alignment(cfg) as reader:
        return reader.header.to_dict()


def list_contigs(reader, cfg) -> list[str]:
    """Resolve which contigs to process from an open reader (config selection)."""
    return cfg.select_contigs(reader.references)


def contig_lengths(reader, contigs) -> dict[str, int]:
    return {c: reader.get_reference_length(c) for c in contigs}


def lengths_from_reference(cfg, contigs) -> dict[str, int]:
    """Contig lengths from the reference .fai alone (no CRAM) -- the CRAM-free
    gather/plot flow. Contigs absent from the reference are dropped."""
    fa = pysam.FastaFile(str(cfg.reference))
    try:
        names = set(fa.references)
        return {c: fa.get_reference_length(c) for c in contigs if c in names}
    finally:
        fa.close()
