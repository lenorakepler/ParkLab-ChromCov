"""
Per-base bedgraph I/O: the write -> read round-trip must reconstruct the depth
vector exactly (it's the resume checkpoint + graph input), including all-zero and
nonzero-at-the-ends cases the RLE encoding is easy to get wrong.
"""
import gzip

import numpy as np

from chromcov.reduce import ChromDepth
from chromcov.io import track as perbase


def _roundtrip(tmp_path, vec):
    perbase.write_bedgraph(tmp_path, "chrT", ChromDepth(vec))
    return perbase.read_bedgraph(perbase.bedgraph_path(tmp_path, "chrT"), vec.size)


def test_roundtrip_exact(tmp_path):
    vec = np.array([0, 0, 3, 3, 3, 0, 7, 7, 0], dtype=np.int32)
    assert np.array_equal(_roundtrip(tmp_path, vec), vec)


def test_roundtrip_all_zero(tmp_path):
    # RLE omits zero runs -> an empty bedgraph -> length still reconstructs as zeros.
    vec = np.zeros(5, dtype=np.int32)
    assert np.array_equal(_roundtrip(tmp_path, vec), vec)


def test_roundtrip_ends_nonzero(tmp_path):
    vec = np.array([4, 4, 1, 1, 1, 9], dtype=np.int32)   # nonzero at both ends
    assert np.array_equal(_roundtrip(tmp_path, vec), vec)


def test_has_and_list(tmp_path):
    assert not perbase.has_bedgraph(tmp_path, "chrT")
    assert perbase.bedgraph_chroms(tmp_path) == []
    perbase.write_bedgraph(tmp_path, "chr1", ChromDepth(np.array([1, 1, 1], dtype=np.int32)))
    perbase.write_bedgraph(tmp_path, "chr2", ChromDepth(np.array([2, 2], dtype=np.int32)))
    assert perbase.has_bedgraph(tmp_path, "chr1")
    assert perbase.bedgraph_chroms(tmp_path) == ["chr1", "chr2"]


def test_write_is_atomic_no_part_left(tmp_path):
    # A successful write leaves the final file and no .part temp behind.
    perbase.write_bedgraph(tmp_path, "chr1", ChromDepth(np.array([5, 5], dtype=np.int32)))
    assert perbase.bedgraph_path(tmp_path, "chr1").exists()
    assert list(tmp_path.glob("*.part")) == []


def test_bedgraph_is_valid_gzip_rle(tmp_path):
    vec = np.array([0, 2, 2, 0, 5], dtype=np.int32)
    perbase.write_bedgraph(tmp_path, "chrT", ChromDepth(vec))
    with gzip.open(perbase.bedgraph_path(tmp_path, "chrT"), "rt") as fh:
        lines = [ln.split("\t") for ln in fh.read().splitlines()]
    # RLE: only the nonzero runs [1,3)=2 and [4,5)=5 are written.
    assert lines == [["chrT", "1", "3", "2"], ["chrT", "4", "5", "5"]]
