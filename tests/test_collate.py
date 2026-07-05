"""
Level 3: comparing analysis runs. Builds a small out/<coverage-key>/<slug>/ tree
(coverage.tsv + run.json) and checks find_runs + pivot + param parsing.
"""
import json

from chromcov.collate import find_runs, pivot


def _make_run(root, cov_key, slug, rows, sidecar):
    d = root / cov_key / slug
    d.mkdir(parents=True)
    cols = list(rows[0].keys())
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
    (d / "coverage.tsv").write_text("\n".join(lines) + "\n")
    (d / "run.json").write_text(json.dumps(sidecar))
    return d


def _sidecar(baseline, strata):
    return {"coverage_key": "covKEY",
            "config": {"coverage": {"min_mapping_quality": 0},
                       "analysis": {"window": 10000, "baseline": baseline, "strata": strata}}}


def test_find_and_pivot_copy_number(tmp_path):
    out = tmp_path / "out"
    _make_run(out, "covKEY", "w10000-nostrata-aaa",
              [{"chrom": "chr1", "length": 100, "bases": 200, "mean_coverage": 2.0, "copy_number": 2.0},
               {"chrom": "chr2", "length": 100, "bases": 300, "mean_coverage": 3.0, "copy_number": 3.0}],
              _sidecar("autosomal-median", {}))
    _make_run(out, "covKEY", "w10000-strata_easy-bbb",
              [{"chrom": "chr1", "length": 100, "bases": 200, "mean_coverage": 2.0, "copy_number": 1.5},
               {"chrom": "chr2", "length": 100, "bases": 300, "mean_coverage": 3.0, "copy_number": 2.5}],
              _sidecar("easy-autosomal-median", {"easy": "e.bed"}))

    runs = find_runs(out)
    assert len(runs) == 2

    run_ids, table = pivot(runs, "copy_number")
    assert len(run_ids) == 2
    assert list(table) == ["chr1", "chr2"]                 # karyotypic order
    assert set(table["chr1"].values()) == {"2.0", "1.5"}   # CN differs across the two runs

    stratified = next(r for r in runs if r.params["strata"] == ["easy"])
    assert stratified.params["baseline"] == "easy-autosomal-median"


def test_pivot_skips_missing_metric(tmp_path):
    out = tmp_path / "out"
    _make_run(out, "k", "w10000-nostrata-ccc",   # a --fast-style run: mean only
              [{"chrom": "chr1", "length": 100, "bases": 200, "mean_coverage": 2.0}],
              {"coverage_key": "k", "config": {"coverage": {}, "analysis": {}}})
    runs = find_runs(out)
    _, table = pivot(runs, "copy_number")
    assert table == {}   # no copy_number column -> no cells
