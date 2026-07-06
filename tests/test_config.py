"""
Config.load: the one place a run is assembled. Config file is the base; the CLI
contributes only overrides. These pin that precedence + the required-inputs rule.
"""
import pytest
import yaml

from chromcov.config import Config, _YAML_TO_FIELD


def test_every_config_field_is_loadable_from_yaml():
    """Guard against the defaults drifting again: every model field must be
    reachable through from_yaml -- either mapped in _YAML_TO_FIELD or one of the
    nested-dict special cases. Add a Config field without wiring it -> this fails."""
    special = {"qc", "strata", "scatter_min_easy_frac"}
    assert set(Config.model_fields) == set(_YAML_TO_FIELD.values()) | special


def _write_cfg(tmp_path, inputs_extra=None, **sections):
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")
    inputs = {"cram": str(cram), "reference": str(ref), **(inputs_extra or {})}
    data = {"inputs": inputs, **sections}
    path = tmp_path / "run.yaml"
    path.write_text(yaml.safe_dump(data))
    return path


def test_cli_override_wins_over_config(tmp_path):
    cfg = _write_cfg(tmp_path, filters={"min_mapping_quality": 5}, analysis={"window": 20_000})
    run = Config.load(cfg, overrides={"min_mapping_quality": 30})
    assert run.min_mapping_quality == 30   # CLI override wins
    assert run.window == 20_000            # untouched value comes from the config


def test_config_only_no_cli(tmp_path):
    cfg = _write_cfg(tmp_path, filters={"min_mapping_quality": 12})
    run = Config.load(cfg)
    assert run.min_mapping_quality == 12
    assert run.window == 10_000            # model default when config is silent


def test_overrides_only_no_config(tmp_path):
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")
    run = Config.load(None, overrides={"cram": str(cram), "reference": str(ref)})
    assert run.min_mapping_quality == 0    # all defaults from the model
    assert run.window == 10_000


def test_requires_inputs():
    with pytest.raises(Exception):   # ValidationError: cram/reference missing
        Config.load(None)


def test_default_index_derived_from_cram(tmp_path):
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")
    run = Config.load(None, overrides={"cram": str(cram), "reference": str(ref)})
    assert run.index == cram.with_suffix(".cram.crai")


def test_chroms_from_config_file(tmp_path):
    cfg = _write_cfg(tmp_path, contigs={"chroms": ["chr1", "chr2"]})
    run = Config.load(cfg)
    assert run.chroms == ("chr1", "chr2")   # explicit subset now settable in YAML
    assert run.select_contigs(["chr1", "chr2", "chr3"]) == ["chr1", "chr2"]


def test_verify_reference_from_config_file(tmp_path):
    cfg = _write_cfg(tmp_path, inputs_extra={"verify_reference": "skip"})
    run = Config.load(cfg)
    assert run.verify_reference == "skip"


def test_scatter_cap_cn_from_config_file(tmp_path):
    cfg = _write_cfg(tmp_path, copy_number={"scatter_cap_cn": 4.0})
    run = Config.load(cfg)
    assert run.scatter_cap_cn == 4.0


def test_qc_thresholds_from_yaml(tmp_path):
    cfg = _write_cfg(tmp_path, qc={"min_median": 8, "gain_cn": 3.0})
    run = Config.load(cfg)
    assert run.qc.min_median == 8
    assert run.qc.gain_cn == 3.0
    assert run.qc.loss_cn == 1.5   # unset key keeps the model default
