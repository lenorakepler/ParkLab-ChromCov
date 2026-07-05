"""
Per-base track store (Level 1) + analysis-run keys (Level 2), no CRAM needed.

Pins the two properties the reuse relies on: the per-base RLE round-trips exactly
(store -> load reconstructs the vector), and the content-addressed keys separate
what should be separate (filter params -> coverage-key; window/strata -> analysis-key).
"""
import numpy as np

from chromcov.analysis import ChromDepth
from chromcov.config import AnalysisConfig, CoverageConfig
from chromcov.perbase import PerBaseStore, analysis_key, analysis_slug


def _cfg(tmp_path, **kw):
    # _stable_input_id stats the files, so they must exist (contents irrelevant).
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"cram-bytes")
    (tmp_path / "x.cram.crai").write_bytes(b"crai")   # provenance stats the index too
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"ref-bytes")
    return CoverageConfig(cram=cram, reference=ref, **kw)


def test_perbase_roundtrip_exact(tmp_path):
    store = PerBaseStore(tmp_path / "perbase", _cfg(tmp_path))
    vec = np.array([0, 0, 3, 3, 3, 0, 7, 7, 0], dtype=np.int32)
    store.store("chr1", ChromDepth(vec))
    assert store.has("chr1")
    assert np.array_equal(store.load("chr1", vec.size), vec)


def test_perbase_roundtrip_all_zero(tmp_path):
    store = PerBaseStore(tmp_path / "perbase", _cfg(tmp_path))
    vec = np.zeros(5, dtype=np.int32)
    store.store("chrZ", ChromDepth(vec))
    assert np.array_equal(store.load("chrZ", vec.size), vec)


def test_perbase_roundtrip_ends_nonzero(tmp_path):
    store = PerBaseStore(tmp_path / "perbase", _cfg(tmp_path))
    vec = np.array([4, 4, 1, 1, 1, 9], dtype=np.int32)  # nonzero at both ends
    store.store("chrE", ChromDepth(vec))
    assert np.array_equal(store.load("chrE", vec.size), vec)


def test_coverage_key_depends_on_filter_params(tmp_path):
    c1 = _cfg(tmp_path)
    c2 = _cfg(tmp_path, min_mapping_quality=20)
    assert PerBaseStore(tmp_path / "pb", c1).key != PerBaseStore(tmp_path / "pb", c2).key


def test_analysis_key_separates_window_and_strata():
    base = analysis_key("COV", AnalysisConfig(), [])
    assert analysis_key("COV", AnalysisConfig(), ["easy", "difficult"]) != base   # strata matters
    assert analysis_key("COV", AnalysisConfig(window=5000), []) != base           # window matters
    assert analysis_key("OTHER", AnalysisConfig(), []) != base                    # coverage-key matters
    # order-independent strata labels -> same key
    assert analysis_key("COV", AnalysisConfig(), ["easy", "difficult"]) == \
           analysis_key("COV", AnalysisConfig(), ["difficult", "easy"])


def test_analysis_slug_is_human_scannable():
    assert "nostrata" in analysis_slug(AnalysisConfig(), [], "abc123")
    slug = analysis_slug(AnalysisConfig(), ["easy", "difficult"], "abc123")
    assert slug.startswith("w10000-strata_difficult_easy") and slug.endswith("abc123")


def test_sidecar_summary_roundtrip(tmp_path):
    store = PerBaseStore(tmp_path / "pb", _cfg(tmp_path))
    store.store("chr1", ChromDepth(np.array([1, 1, 1], dtype=np.int32)))
    store.write_sidecar({"chr1": {"length": 3, "bases": 3, "mean": 1.0}})
    assert store.exists()
    assert store.chroms() == ["chr1"]
    assert store.read_summary()["chr1"]["mean"] == 1.0
