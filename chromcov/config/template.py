"""
Generate an editable run-config YAML from the *live* Config model.

The template is never hand-maintained: section grouping, per-field comments, and
every value are read from `Config` / `QCThresholds` model fields and their
`json_schema_extra` metadata at generation time (was gen_config.py). Adding a
Config field with its `section`/`comment` metadata makes it appear here for free,
so an emitted config always matches the current code defaults.

`chromcov gen-config -o my.yaml` writes it; the layout mirrors what
`Config.from_yaml` reads back, so the file round-trips.
"""
from __future__ import annotations

from pathlib import Path

from ..filtering import SAM_FLAGS
from .schema import Config, QCThresholds
from .schema import _meta

# Regular sections rendered by walking the model metadata. The three irregular
# fields (qc / strata / scatter_min_easy_frac) are rendered by hand below.
_REGULAR_SECTIONS = ["inputs", "filters", "contigs", "analysis", "copy_number", "output"]

_SECTION_HEADINGS = {
    "inputs": "inputs -- default to the bundled Park Lab COLO829T test data",
    "filters": "read filtering (samtools flag/MAPQ semantics)",
    "contigs": "contig selection -- include globs first, then exclude",
    "analysis": "analysis (used only with --full)",
    "copy_number": "copy number (used only with --full)",
    "output": "output",
}


def _config_defaults() -> dict:
    """Resolved default value per Config field (calls default_factory; required
    fields -> None so they render as commented placeholders)."""
    out: dict = {}
    for name, fi in Config.model_fields.items():
        if fi.is_required():
            out[name] = None
        elif fi.default_factory is not None:
            out[name] = fi.default_factory()  # type: ignore[call-arg]
        else:
            out[name] = fi.default
    return out


def _comment(text: str, indent: str = "  ") -> str:
    return "\n".join(f"{indent}# {ln}" for ln in text.split("\n"))


def _inline(seq) -> str:
    """A sequence as an inline YAML list, quoting strings, bare numbers."""
    parts = [f'"{x}"' if isinstance(x, str) else str(x) for x in seq]
    return "[" + ", ".join(parts) + "]"


def _scalar(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return _inline(value)
    if isinstance(value, Path):
        return str(value)
    return str(value)


def _flag_names(mask: int) -> list[str]:
    """Bitmask -> the SAM flag names that compose it (config-friendly form)."""
    return [name for name, bit in SAM_FLAGS.items() if mask & bit]


def _flag_value(key: str, mask: int) -> str:
    """A -f/-F/-G flags entry: inline `[]` when empty, else a YAML block list."""
    names = _flag_names(mask)
    if not names:
        return f"  {key}: []"
    return "\n".join([f"  {key}:"] + [f"    - {n}" for n in names])


def _render_field(name, fi, value) -> str:
    m = _meta(fi)
    key = m.get("yaml_key", name)
    comment = m.get("comment", "")
    lines = [_comment(comment)] if comment else []

    if m.get("kind") == "flags":
        lines.append(_flag_value(key, value))
    elif fi.is_required():
        lines.append(f"  # {key}: {m.get('placeholder', '...')}")
    elif value is None:
        if "placeholder" in m:
            lines.append(f"  # {key}: {m['placeholder']}")
        else:
            lines.append(f"  {key}: null")
    else:
        lines.append(f"  {key}: {_scalar(value)}")
    return "\n".join(lines)


def _render_section(section: str, defaults: dict) -> str:
    fields = [(n, fi) for n, fi in Config.model_fields.items()
              if _meta(fi).get("section") == section]
    heading = _SECTION_HEADINGS.get(section, section)
    body = "\n".join(_render_field(n, fi, defaults[n]) for n, fi in fields)
    return f"# --- {heading} ---\n{section}:\n{body}"


def _render_qc() -> str:
    qc = QCThresholds()
    lines = ["# --- QC abnormality-flag thresholds (set any subset; unset keys default) ---", "qc:"]
    for name, fi in QCThresholds.model_fields.items():
        comment = _meta(fi).get("comment", "")
        val = str(getattr(qc, name))
        lines.append(f"  {name}: {val:<8}# {comment}" if comment else f"  {name}: {val}")
    return "\n".join(lines)


def _render_strata(defaults: dict) -> str:
    return (
        "# --- callability strata (SMaHT_Regional_Categorization, GRCh38) ---\n"
        "# label -> BED[.gz] path. Empty by default. Fetch with `chromcov fetch strata`,\n"
        "# then fill these in or use --strata.\n"
        "strata: {}\n"
        "  # easy: data/SMaHT_easy_hg38.bed.gz\n"
        "  # difficult: data/SMaHT_difficult_hg38.bed.gz\n"
        "  # extreme: data/SMaHT_extreme_hg38.bed.gz\n"
        f"  # scatter_min_easy_frac: {defaults['scatter_min_easy_frac']}   "
        "# >0 drops windows mostly outside 'easy' from the scatter"
    )


def render_default_config() -> str:
    """Build the commented config text from the live model defaults + metadata."""
    d = _config_defaults()
    header = (
        "# chromcov run configuration -- every option with its current default value.\n"
        "#\n"
        "# GENERATED by `chromcov gen-config` from the live Config model, so it always\n"
        "# matches the code. Inputs default to the bundled Park Lab COLO829T test data\n"
        "# (get it with `chromcov fetch inputs`); point inputs.cram / inputs.reference\n"
        "# elsewhere to run on your own. CLI flags override anything here.\n"
        f"# Valid flag names: {', '.join(SAM_FLAGS)}. Integer bitmasks also accepted."
    )
    blocks = [header]
    blocks += [_render_section(s, d) for s in ["inputs", "filters", "contigs", "analysis", "copy_number"]]
    blocks.append(_render_qc())
    blocks.append(_render_strata(d))
    blocks.append(_render_section("output", d))
    return "\n\n".join(blocks) + "\n"


def write_default_config(path: str | Path) -> Path:
    """Render the live-default config and write it to `path`."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(render_default_config())
    return path
