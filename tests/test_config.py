"""
Config.load: the one place a run is assembled. Config file is the base; the CLI
contributes only overrides. These pin that precedence + the required-inputs rule.
"""
import pytest
import yaml

from chromcov.config import Config


def _write_cfg(tmp_path, **sections):
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")
    data = {"inputs": {"cram": str(cram), "reference": str(ref)}, **sections}
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


def test_qc_thresholds_from_yaml(tmp_path):
    cfg = _write_cfg(tmp_path, qc={"min_median": 8, "gain_cn": 3.0})
    run = Config.load(cfg)
    assert run.qc.min_median == 8
    assert run.qc.gain_cn == 3.0
    assert run.qc.loss_cn == 1.5   # unset key keeps the model default
