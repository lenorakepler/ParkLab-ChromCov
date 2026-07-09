"""
Run configuration -- the SOLE source of field truth.

Every option lives here once, as a Pydantic field, together with the metadata
(`section` / `yaml_key` / `comment`, in `json_schema_extra`) that both the YAML
template generator (`config/template.py`) and the flat YAML<->field map
(`_YAML_TO_FIELD`, derived below) read. Add a field with its section metadata and
it shows up in the generated config and the loader for free -- nothing is
hand-maintained in two places.

Read filtering (SAM flags, `ReadFilter`, `to_mask`) lives in `chromcov.filtering`.
"""
from __future__ import annotations

import fnmatch
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from ..filtering import (
    DEFAULT_EXCLUDE,
    DEFAULT_EXCLUDE_CONTIGS,
    DEFAULT_INCLUDE_CONTIGS,
    to_mask,
)


class QCThresholds(BaseModel):
    """User-tunable thresholds for the coverage-QC abnormality flags (the logic is
    in policy.py). Set any subset in a --config YAML under a `qc:` section, e.g.
    `qc: {min_median: 8, min_breadth_20x: 0.5}`; unset keys keep these defaults."""

    model_config = ConfigDict(extra="forbid")

    gain_cn: float = Field(2.5, json_schema_extra={"comment": "CN >= this (autosomes) -> CN_GAIN"})
    loss_cn: float = Field(1.5, json_schema_extra={"comment": "CN <= this (autosomes) -> CN_LOSS"})
    depleted_cn: float = Field(0.25, json_schema_extra={"comment": "CN <= this (any contig) -> CN_DEPLETED"})
    min_median: float = Field(10.0, json_schema_extra={"comment": "median depth < this -> LOW_DEPTH"})
    uneven_robust_cv: float = Field(0.5, json_schema_extra={"comment": "MAD/median > this -> UNEVEN"})
    min_breadth_20x: float = Field(0.70, json_schema_extra={"comment": "breadth@20x < this -> LOW_CALLABLE"})
    extreme_median_mult: float = Field(5.0, json_schema_extra={"comment": "median > this * baseline -> EXTREME_DEPTH"})


def _meta(fi) -> dict:
    """The json_schema_extra metadata dict for a field (empty if none)."""
    extra = fi.json_schema_extra
    return extra if isinstance(extra, dict) else {}


class Config(BaseModel):
    """A whole run: coverage knobs + the extension knobs `--full` adds. Validated
    at construction (the CLI/YAML boundary).

    Each field's `json_schema_extra` carries its YAML `section`, optional
    `yaml_key` (defaults to the field name), and a `comment`. `_YAML_TO_FIELD` and
    the `gen-config` template are both derived from this metadata."""

    model_config = ConfigDict(extra="forbid")

    # --- inputs ---
    cram: Path = Field(..., json_schema_extra={
        "section": "inputs", "placeholder": "path/to/sample.cram",
        "comment": "REQUIRED -- the CRAM to measure"})
    reference: Path = Field(..., json_schema_extra={
        "section": "inputs", "placeholder": "path/to/genome.fa",
        "comment": "REQUIRED -- FASTA the CRAM was compressed against"})
    index: Path | None = Field(None, json_schema_extra={
        "section": "inputs", "comment": "CRAM index; null -> <cram>.crai"})
    verify_reference: Literal["auto", "full", "skip"] = Field("auto", json_schema_extra={
        "section": "inputs", "comment": "reference-M5 preflight check: auto | full | skip"})

    # --- contig selection ---
    # Optional explicit subset (None = all); else include-then-exclude globs over
    # the reference names, to keep decoy/unplaced artifact contigs out of the report.
    chroms: tuple[str, ...] | None = Field(None, json_schema_extra={
        "section": "contigs", "comment": "explicit subset; overrides include/exclude (default: all)",
        "placeholder": "[chr1, chr2]"})
    include_contigs: tuple[str, ...] = Field(DEFAULT_INCLUDE_CONTIGS, json_schema_extra={
        "section": "contigs", "yaml_key": "include", "comment": "keep contigs matching these globs"})
    exclude_contigs: tuple[str, ...] = Field(DEFAULT_EXCLUDE_CONTIGS, json_schema_extra={
        "section": "contigs", "yaml_key": "exclude", "comment": "then drop contigs matching these"})

    # --- read filtering (samtools semantics: -Q/-f/-F/-G) ---
    min_mapping_quality: int = Field(0, json_schema_extra={
        "section": "filters", "comment": "-Q: minimum mapping quality a read must have"})
    include_flags: int = Field(0, json_schema_extra={
        "section": "filters", "kind": "flags", "comment": "-f: require ALL of these flags"})
    exclude_flags: int = Field(DEFAULT_EXCLUDE, json_schema_extra={
        "section": "filters", "kind": "flags",
        "comment": "-F: exclude a read if it has ANY of these flags\n"
                   "`duplicate` is deliberately NOT excluded by default -- in high-depth cancer\n"
                   "data, coordinate collisions make duplicate-marking drop real coverage; add\n"
                   "`- duplicate` here to reproduce a standard dedup'd number."})
    exclude_all_flags: int = Field(0, json_schema_extra={
        "section": "filters", "kind": "flags", "comment": "-G: exclude a read ONLY if it has ALL of these flags"})

    # --- extension knobs (ignored unless --full) ---
    window: int = Field(10_000, json_schema_extra={
        "section": "analysis", "comment": "windowed-mean bin size (bp)"})
    hist_cap: int = Field(200_000, json_schema_extra={
        "section": "analysis", "comment": "clip extreme pileup depths into a top bin"})
    breadth_thresholds: tuple[int, ...] = Field((1, 5, 10, 15, 20, 30), json_schema_extra={
        "section": "analysis", "comment": "depth cutoffs for breadth-of-coverage"})
    ploidy: int = Field(2, json_schema_extra={"section": "copy_number", "comment": "assumed baseline ploidy"})
    baseline: Literal["easy-autosomal-median", "autosomal-median"] = Field(
        "easy-autosomal-median", json_schema_extra={
            "section": "copy_number", "comment": "easy-autosomal-median | autosomal-median"})
    scatter_cap_cn: float = Field(6.0, json_schema_extra={
        "section": "copy_number", "comment": "cap on approx copy number on the scatter y-axis"})

    # --- irregular fields (special-cased in from_yaml + the template) ---
    qc: QCThresholds = Field(default_factory=QCThresholds)   # abnormality-flag thresholds
    strata: dict[str, str] = {}          # label -> BED[.gz] path
    # 0 -> scatter shows all windows colored by callability tier; >0 -> restrict
    # to callable ('easy') windows (the old callable-only view).
    scatter_min_easy_frac: float = 0.0

    # --- output ---
    outdir: Path = Field(Path("out"), json_schema_extra={"section": "output", "comment": "output root"})
    plots: bool = Field(True, json_schema_extra={"section": "output", "comment": "render plots (--full)"})

    @field_validator("include_flags", "exclude_flags", "exclude_all_flags", mode="before")
    @classmethod
    def _normalize_mask(cls, v):
        return to_mask(v)

    @model_validator(mode="after")
    def _default_index(self) -> "Config":
        if self.index is None and self.cram is not None:
            self.index = self.cram.with_suffix(self.cram.suffix + ".crai")
        return self

    def select_contigs(self, references) -> list[str]:
        """Resolve which contigs to process: explicit `chroms` wins, else apply
        include-then-exclude globs over the CRAM's reference names."""
        if self.chroms is not None:
            return list(self.chroms)
        kept = [r for r in references
                if any(fnmatch.fnmatch(r, p) for p in self.include_contigs)]
        return [r for r in kept
                if not any(fnmatch.fnmatch(r, p) for p in self.exclude_contigs)]

    @classmethod
    def from_yaml(cls, path: str | Path) -> dict:
        """Flatten a run-config YAML (see gen-config output) into a dict of only the
        `Config` field names the file actually sets -- absent keys are omitted so the
        model supplies the default (defaults live in one place: the model). Returned
        as a dict (not a `Config`) so `load` can layer CLI overrides on top before
        validating once. Type coercion (list->tuple, flag names->mask) is the model's
        job, not this function's."""
        data = yaml.safe_load(Path(path).read_text()) or {}

        payload: dict = {}
        for (section, key), field in _YAML_TO_FIELD.items():
            sec = data.get(section) or {}
            if key in sec and sec[key] is not None:
                payload[field] = sec[key]

        # qc: a partial dict is fine -- QCThresholds fills the unset keys.
        if data.get("qc"):
            payload["qc"] = data["qc"]

        # strata: the {label: path} map plus the scatter knob that lives beside it.
        strata = dict(data.get("strata") or {})
        if "scatter_min_easy_frac" in strata:
            payload["scatter_min_easy_frac"] = strata.pop("scatter_min_easy_frac")
        if strata:
            payload["strata"] = strata

        return payload

    @classmethod
    def load(cls, path: str | Path | None = None, overrides: dict | None = None) -> "Config":
        """Assemble a run: the config file (if any) is the base, the CLI passes a
        dict of only the options it was explicitly given, merged on top. Defaults
        and validation live in the model -- nothing is duplicated in the CLI. This
        is the whole 'config authoritative, CLI overrides' rule, in one place."""
        base = cls.from_yaml(path) if path else {}
        return cls.model_validate({**base, **(overrides or {})})


# YAML (section, key) -> flat Config field, DERIVED from the field metadata above
# (no hand-maintained table). Every field carrying a `section` is reachable through
# from_yaml; the three irregular fields (qc / strata / scatter_min_easy_frac) have
# no section and are handled explicitly in from_yaml.
_YAML_TO_FIELD = {
    (_meta(fi)["section"], _meta(fi).get("yaml_key", name)): name
    for name, fi in Config.model_fields.items()
    if "section" in _meta(fi)
}
