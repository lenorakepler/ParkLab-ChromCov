"""
Shared coverage configuration (Pydantic v2).

The whole point of the dispatcher is that BOTH backends (the hand-rolled pysam
calculator and mosdepth) consume this one object. Anything backend-specific
lives behind the backend boundary, not here. Where a knob can't be honored by a
backend, that backend should raise -- silent divergence between backends is the
worst outcome, since the main reason to have both is cross-validation.

Config is the trust boundary: `CoverageConfig` and `AnalysisConfig` are Pydantic
models so a run driven by a YAML file or CLI args is *validated* (str->Path
coercion, flag-name-list -> int-mask normalization, unknown-backend rejection)
with actionable errors, instead of failing deep in the fetch loop. The computed
value objects (`ChromCoverage`, `ChromStats`) stay plain dataclasses -- they're
internal results, not parsed input.

Flag masks mirror samtools/mosdepth integer semantics. Note the default-mask
mismatch between tools, spelled out below: we pin an explicit shared default so
the two backends actually agree out of the box.
"""
from __future__ import annotations

import fnmatch
from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

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

# Explicit shared default = unmapped | secondary | qcfail | duplicate | supplementary.
# This is deliberately the UNION of the two tools' native defaults so neither
# backend is silently doing something different:
#   - read_filter.py default   = unmapped|secondary|duplicate|supplementary  (3332)
#   - mosdepth --flag default   = unmapped|secondary|qcfail|duplicate        (1796)
# Pinning it here (3844) removes that footgun; override per-run if you want to
# reproduce a specific tool's out-of-the-box number.
DEFAULT_EXCLUDE = (
    SAM_FLAGS["unmapped"]
    | SAM_FLAGS["secondary"]
    | SAM_FLAGS["qcfail"]
    | SAM_FLAGS["duplicate"]
    | SAM_FLAGS["supplementary"]
)

DEFAULT_INCLUDE_CONTIGS = ("chr*",)
DEFAULT_EXCLUDE_CONTIGS = ("*_alt", "*_random", "chrUn*", "*_decoy", "HLA*", "chrEBV")


def to_mask(flags) -> int:
    """Normalize a flag spec (int mask | iterable of names/ints | None) to an int.

    Raises KeyError-free ValueError on an unknown flag name so a typo in a YAML
    config surfaces as a clear message rather than a bare KeyError.
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


class CoverageConfig(BaseModel):
    """The knobs both backends share. Validated at construction (CLI/YAML boundary)."""

    model_config = ConfigDict(extra="forbid")

    # --- inputs (shared) ---
    cram: Path
    reference: Path
    index: Path | None = None            # defaults to <cram>.crai if None

    # --- backend selection ---
    backend: Literal["native", "mosdepth"] = "native"

    # Optional explicit contig subset (None = all). native loops only these;
    # mosdepth would need a --chrom (single) or a BED to match.
    chroms: tuple[str, ...] | None = None

    # Glob include/exclude over reference names, applied when `chroms` is None.
    # Keeps decoy/unplaced artifact contigs out of the headline report.
    include_contigs: tuple[str, ...] = DEFAULT_INCLUDE_CONTIGS
    exclude_contigs: tuple[str, ...] = DEFAULT_EXCLUDE_CONTIGS

    # --- read filtering (shared knobs; samtools/mosdepth semantics) ---
    min_mapping_quality: int = 0         # native -Q  /  mosdepth --mapq
    include_flags: int = 0               # native -f  /  mosdepth --include-flag
    exclude_flags: int = DEFAULT_EXCLUDE  # native -F  /  mosdepth --flag
    exclude_all_flags: int = 0           # native -G  /  mosdepth: NO EQUIVALENT (backend raises)

    # --- outputs / execution ---
    per_base: bool = False               # keep full per-base depth vector (native) /
                                         #   drop mosdepth --no-per-base when True
    threads: int = 1                     # mosdepth --threads; native is single-threaded

    @field_validator("include_flags", "exclude_flags", "exclude_all_flags", mode="before")
    @classmethod
    def _normalize_mask(cls, v):
        # Accept an int mask, a list of flag names/ints, or None -> int mask.
        return to_mask(v)

    @model_validator(mode="after")
    def _default_index(self) -> "CoverageConfig":
        if self.index is None:
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
    def from_yaml(cls, path: str | Path) -> "CoverageConfig":
        """Build a CoverageConfig from the inputs/backend/filters/contigs sections
        of a run-config YAML (see config.example.yml). Pipeline-level sections
        (analysis/copy_number/strata/output) are consumed by AnalysisConfig."""
        data = yaml.safe_load(Path(path).read_text()) or {}
        inputs = data.get("inputs", {}) or {}
        backend = data.get("backend", {}) or {}
        filters = data.get("filters", {}) or {}
        contigs = data.get("contigs", {}) or {}
        payload = {
            "cram": inputs.get("cram"),
            "reference": inputs.get("reference"),
            "index": inputs.get("index"),
            "backend": backend.get("name", "native"),
            "threads": backend.get("threads", 1),
            "min_mapping_quality": filters.get("min_mapping_quality", 0),
            "include_flags": filters.get("include_flags", 0),
            "exclude_flags": filters.get("exclude_flags", DEFAULT_EXCLUDE),
            "exclude_all_flags": filters.get("exclude_all_flags", 0),
        }
        if "include" in contigs:
            payload["include_contigs"] = tuple(contigs["include"])
        if "exclude" in contigs:
            payload["exclude_contigs"] = tuple(contigs["exclude"])
        return cls.model_validate(payload)


class AnalysisConfig(BaseModel):
    """Pipeline-level knobs: everything downstream of the per-base pass. Kept
    separate from CoverageConfig because these don't affect the coverage numbers,
    only how they're reduced/reported."""

    model_config = ConfigDict(extra="forbid")

    # --- windowed track / histogram ---
    window: int = 10_000
    hist_cap: int = 200_000
    breadth_thresholds: tuple[int, ...] = (1, 5, 10, 15, 20, 30)

    # --- copy number ---
    ploidy: int = 2
    baseline: Literal["easy-autosomal-median", "autosomal-median"] = "easy-autosomal-median"

    # --- callability strata: label -> BED[.gz] path ---
    strata: dict[str, str] = {}
    scatter_min_easy_frac: float = 0.5

    # --- output ---
    outdir: Path = Path("out")
    per_base: bool = False               # write the RLE bedgraph (opt-in; bulky on WGS)
    plots: bool = True

    @classmethod
    def from_yaml(cls, path: str | Path) -> "AnalysisConfig":
        data = yaml.safe_load(Path(path).read_text()) or {}
        analysis = data.get("analysis", {}) or {}
        copy_number = data.get("copy_number", {}) or {}
        strata = dict(data.get("strata", {}) or {})
        output = data.get("output", {}) or {}
        scatter_min = strata.pop("scatter_min_easy_frac", 0.5)
        per_base = (output.get("per_base", {}) or {}).get("enabled", False)
        payload = {
            "window": analysis.get("window", 10_000),
            "hist_cap": analysis.get("hist_cap", 200_000),
            "breadth_thresholds": tuple(analysis.get("breadth_thresholds", (1, 5, 10, 15, 20, 30))),
            "ploidy": copy_number.get("ploidy", 2),
            "baseline": copy_number.get("baseline", "easy-autosomal-median"),
            "strata": strata,
            "scatter_min_easy_frac": scatter_min,
            "outdir": output.get("outdir", "out"),
            "per_base": per_base,
            "plots": output.get("plots", True),
        }
        return cls.model_validate(payload)


class RunConfig(BaseModel):
    """A whole run = coverage config + analysis config.

    The single place a run is assembled from, so the "config is authoritative,
    CLI only overrides" rule lives in one method (`load`) instead of being
    re-implemented per command. The config file (if given) supplies the base; the
    CLI passes dicts of only the options it was *explicitly* given, which are
    merged on top. Defaults and validation stay in the two models -- nothing is
    duplicated in the CLI."""

    model_config = ConfigDict(extra="forbid")

    coverage: CoverageConfig
    analysis: AnalysisConfig = AnalysisConfig()

    @classmethod
    def load(cls, path: str | Path | None = None, *,
             coverage: dict | None = None, analysis: dict | None = None) -> "RunConfig":
        cov_base = CoverageConfig.from_yaml(path).model_dump() if path else {}
        ana_base = AnalysisConfig.from_yaml(path).model_dump() if path else {}
        return cls(
            coverage=CoverageConfig.model_validate({**cov_base, **(coverage or {})}),
            analysis=AnalysisConfig.model_validate({**ana_base, **(analysis or {})}),
        )
