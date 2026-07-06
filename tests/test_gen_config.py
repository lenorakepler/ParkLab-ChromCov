"""
gen-config must produce a file that (a) parses, and (b) round-trips back to the
exact model defaults -- that's the whole promise ("always matches the code").
If a Config default changes, this test guarantees the generated file follows.
"""
import yaml

from chromcov.config import Config, QCThresholds
from chromcov.gen_config import render_default_config, write_default_config


def test_generated_config_is_valid_yaml():
    data = yaml.safe_load(render_default_config())
    assert set(data) >= {"inputs", "filters", "contigs", "analysis",
                         "copy_number", "qc", "strata", "output"}


def test_roundtrip_matches_model_defaults(tmp_path):
    path = write_default_config(tmp_path / "gen.yaml")
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")

    run = Config.load(path, overrides={"cram": str(cram), "reference": str(ref)})

    defaults = Config.model_fields
    assert run.min_mapping_quality == defaults["min_mapping_quality"].default
    assert run.include_flags == defaults["include_flags"].default
    assert run.exclude_flags == defaults["exclude_flags"].default   # names -> mask
    assert run.exclude_all_flags == defaults["exclude_all_flags"].default
    assert run.include_contigs == defaults["include_contigs"].default
    assert run.exclude_contigs == defaults["exclude_contigs"].default
    assert run.window == defaults["window"].default
    assert run.hist_cap == defaults["hist_cap"].default
    assert run.breadth_thresholds == defaults["breadth_thresholds"].default
    assert run.ploidy == defaults["ploidy"].default
    assert run.baseline == defaults["baseline"].default
    assert run.scatter_min_easy_frac == defaults["scatter_min_easy_frac"].default
    assert run.scatter_cap_cn == defaults["scatter_cap_cn"].default
    assert run.verify_reference == defaults["verify_reference"].default
    assert run.chroms == defaults["chroms"].default   # None -> all contigs
    assert run.plots == defaults["plots"].default
    assert run.strata == {}
    assert run.qc.model_dump() == QCThresholds().model_dump()


def test_tracks_edited_default(tmp_path, monkeypatch):
    """A changed model default flows into the generated file with no template edit."""
    monkeypatch.setattr(Config.model_fields["window"], "default", 25_000)
    data = yaml.safe_load(render_default_config())
    assert data["analysis"]["window"] == 25_000
