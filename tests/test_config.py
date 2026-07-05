"""
RunConfig: the one place a run is assembled. Config file is the base; the CLI
contributes only overrides. These pin that precedence + the required-inputs rule.
"""
import pytest
import yaml

from chromcov.config import RunConfig


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
    run = RunConfig.load(cfg, coverage={"min_mapping_quality": 30})
    assert run.coverage.min_mapping_quality == 30   # CLI override wins
    assert run.analysis.window == 20_000            # untouched value comes from the config


def test_config_only_no_cli(tmp_path):
    cfg = _write_cfg(tmp_path, filters={"min_mapping_quality": 12})
    run = RunConfig.load(cfg)
    assert run.coverage.min_mapping_quality == 12
    assert run.analysis.window == 10_000            # model default when config is silent


def test_overrides_only_no_config(tmp_path):
    cram = tmp_path / "x.cram"
    cram.write_bytes(b"c")
    ref = tmp_path / "r.fa"
    ref.write_bytes(b"r")
    run = RunConfig.load(None, coverage={"cram": str(cram), "reference": str(ref)})
    assert run.coverage.min_mapping_quality == 0    # all defaults from the model
    assert run.analysis.window == 10_000


def test_requires_inputs():
    with pytest.raises(Exception):   # ValidationError: cram/reference missing
        RunConfig.load(None)
