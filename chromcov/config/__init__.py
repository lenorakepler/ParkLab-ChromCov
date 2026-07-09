"""Configuration subpackage: the schema (field-truth source) + YAML template."""
from .schema import Config, QCThresholds, _YAML_TO_FIELD

__all__ = ["Config", "QCThresholds", "_YAML_TO_FIELD"]
