"""
Run Configuration + Read Filtering
"""
from __future__ import annotations

import fnmatch
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

# SAM flag name -> bit (kept here so config/CLI can speak names, internals use ints).
SAM_FLAGS = {
    "paired":        0x1,
    "proper_pair":   0x2,
    "unmapped":      0x4,
    "mate_unmapped": 0x8,
    "reverse":       0x10,
    "mate_reverse":  0x20,
    "read1":         0x40,
    "read2":         0x80,
    "secondary":     0x100,
    "qcfail":        0x200,
    "duplicate":     0x400,
    "supplementary": 0x800,
}

# Default read exclusions: unmapped | secondary | qcfail | supplementary  (2820).
# Duplicates are deliberately KEPT. In high-depth cancer data, duplicate-marking
# flags reads that share start/end coordinates -- but at high depth independent
# fragments increasingly collide on coordinates by chance, so removing "duplicates"
# can under-count genuinely amplified regions. Override with -F to reproduce a
# tool's native default (e.g. mosdepth also excludes duplicates). The mosdepth
# cross-check add-on defaults to THIS same mask, so the two still agree out of the box.
DEFAULT_EXCLUDE = (
    SAM_FLAGS["unmapped"]
    | SAM_FLAGS["secondary"]
    | SAM_FLAGS["qcfail"]
    | SAM_FLAGS["supplementary"]
)

DEFAULT_INCLUDE_CONTIGS = ("chr*",)
DEFAULT_EXCLUDE_CONTIGS = ("*_alt", "*_random", "chrUn*", "*_decoy", "HLA*", "chrEBV")

def to_mask(flags) -> int:
    """
    Outputs a bitmask: either integers, flags (will be converted), or None.
    Quicker filter while still allowing for easy specification in config.

    Raises ValueError on unknown flag
    """
    if flags is None:
        return 0
    if isinstance(flags, int):
        return flags
    mask = 0
    for f in flags:
        if isinstance(f, int):
            mask |= f
        elif f in SAM_FLAGS:
            mask |= SAM_FLAGS[f]
        else:
            raise ValueError(
                f"unknown SAM flag {f!r}; choose from {sorted(SAM_FLAGS)} or pass an int mask"
            )
    return mask

@dataclass
class ReadFilter:
    """
    Build a read filter based on config specification. Create once, use as filter for each read.
    """
    # -f: require a read to have ALL these flags
    include_flags: int = 0  

    # -F: exclude a read if it has ANY of these flags
    exclude_flags: int = DEFAULT_EXCLUDE

    # -G: exclude a read only if it has ALL these flags
    exclude_all_flags: int = 0

    # -Q: minimum mapping quality a read must have
    min_mapping_quality: int = 0

    # Create the mask
    def __post_init__(self):
        self.include_flags     = to_mask(self.include_flags)
        self.exclude_flags     = to_mask(self.exclude_flags)
        self.exclude_all_flags = to_mask(self.exclude_all_flags)

    def fails(self, read):
        flag = read.flag
        if flag & self.exclude_flags:
            return True
        if (flag & self.include_flags) != self.include_flags:
            return True
        if self.exclude_all_flags and (flag & self.exclude_all_flags) == self.exclude_all_flags:
            return True
        if read.mapping_quality < self.min_mapping_quality:
            return True
        return False


class QCThresholds(BaseModel):
    """User-tunable thresholds for the coverage-QC abnormality flags (the logic is
    in qc_flags.py). Set any subset in a --config YAML under a `qc:` section, e.g.
    `qc: {min_median: 8, min_breadth_20x: 0.5}`; unset keys keep these defaults."""

    model_config = ConfigDict(extra="forbid")

    gain_cn: float = 2.5              # CN >= this (autosomes) -> CN_GAIN
    loss_cn: float = 1.5             # CN <= this (autosomes) -> CN_LOSS
    depleted_cn: float = 0.25        # CN <= this (any contig) -> CN_DEPLETED
    min_median: float = 10.0         # median depth < this -> LOW_DEPTH
    uneven_robust_cv: float = 0.5    # MAD/median > this -> UNEVEN
    min_breadth_20x: float = 0.70    # breadth@20x < this -> LOW_CALLABLE
    extreme_median_mult: float = 5.0  # median > this * baseline -> EXTREME_DEPTH


class Config(BaseModel):
    """A whole run: coverage knobs + the extension knobs `--full` adds. Validated
    at construction (the CLI/YAML boundary)."""

    model_config = ConfigDict(extra="forbid")

    # --- inputs ---
    cram: Path
    reference: Path
    index: Path | None = None            # defaults to <cram>.crai if None

    # --- contig selection ---
    # Optional explicit subset (None = all); else include-then-exclude globs over
    # the reference names, to keep decoy/unplaced artifact contigs out of the report.
    chroms: tuple[str, ...] | None = None
    include_contigs: tuple[str, ...] = DEFAULT_INCLUDE_CONTIGS
    exclude_contigs: tuple[str, ...] = DEFAULT_EXCLUDE_CONTIGS

    # --- read filtering (samtools semantics: -Q/-f/-F/-G) ---
    min_mapping_quality: int = 0
    include_flags: int = 0               # -f: require ALL these bits
    exclude_flags: int = DEFAULT_EXCLUDE  # -F: exclude if ANY set
    exclude_all_flags: int = 0           # -G: exclude only if ALL set

    # --- extension knobs (ignored unless --full) ---
    window: int = 10_000
    hist_cap: int = 200_000
    breadth_thresholds: tuple[int, ...] = (1, 5, 10, 15, 20, 30)
    ploidy: int = 2
    baseline: Literal["easy-autosomal-median", "autosomal-median"] = "easy-autosomal-median"
    qc: QCThresholds = Field(default_factory=QCThresholds)   # abnormality-flag thresholds
    strata: dict[str, str] = {}          # label -> BED[.gz] path
    # 0 -> scatter shows all windows colored by callability tier; >0 -> restrict
    # to callable ('easy') windows (the old callable-only view).
    scatter_min_easy_frac: float = 0.0

    # --- output ---
    outdir: Path = Path("out")
    plots: bool = True

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
        """Flatten a run-config YAML (see config.example.yaml) into a dict of
        `Config` field names. Returned as a dict (not a `Config`) so `load` can
        layer CLI overrides on top before validating once."""
        data = yaml.safe_load(Path(path).read_text()) or {}
        inputs = data.get("inputs", {}) or {}
        filters = data.get("filters", {}) or {}
        contigs = data.get("contigs", {}) or {}
        analysis = data.get("analysis", {}) or {}
        copy_number = data.get("copy_number", {}) or {}
        strata = dict(data.get("strata", {}) or {})
        output = data.get("output", {}) or {}
        scatter_min = strata.pop("scatter_min_easy_frac", 0.0)

        payload: dict = {
            "cram": inputs.get("cram"),
            "reference": inputs.get("reference"),
            "index": inputs.get("index"),
            "min_mapping_quality": filters.get("min_mapping_quality", 0),
            "include_flags": filters.get("include_flags", 0),
            "exclude_flags": filters.get("exclude_flags", DEFAULT_EXCLUDE),
            "exclude_all_flags": filters.get("exclude_all_flags", 0),
            "window": analysis.get("window", 10_000),
            "hist_cap": analysis.get("hist_cap", 200_000),
            "breadth_thresholds": tuple(analysis.get("breadth_thresholds", (1, 5, 10, 15, 20, 30, 100, 200))),
            "ploidy": copy_number.get("ploidy", 2),
            "baseline": copy_number.get("baseline", "easy-autosomal-median"),
            "qc": data.get("qc", {}) or {},     # partial dict -> QCThresholds (unset keys default)
            "strata": strata,
            "scatter_min_easy_frac": scatter_min,
            "outdir": output.get("outdir", "out"),
            "plots": output.get("plots", True),
        }
        if "include" in contigs:
            payload["include_contigs"] = tuple(contigs["include"])
        if "exclude" in contigs:
            payload["exclude_contigs"] = tuple(contigs["exclude"])
        # Drop keys the YAML didn't set, so model defaults apply.
        return {k: v for k, v in payload.items() if v is not None}

    @classmethod
    def load(cls, path: str | Path | None = None, overrides: dict | None = None) -> "Config":
        """Assemble a run: the config file (if any) is the base, the CLI passes a
        dict of only the options it was explicitly given, merged on top. Defaults
        and validation live in the model -- nothing is duplicated in the CLI. This
        is the whole 'config authoritative, CLI overrides' rule, in one place."""
        base = cls.from_yaml(path) if path else {}
        return cls.model_validate({**base, **(overrides or {})})
