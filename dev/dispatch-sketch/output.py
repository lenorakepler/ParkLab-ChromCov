"""
Output handling: per-run archival + cross-run comparison.

Two jobs, two shapes (see the design note in the sketch README):

  write_run(...)  -> archival. One directory per run, named by a deterministic
                     hash of (params + input identity), holding a self-describing
                     TSV (# comment header) and a provenance.json sidecar. Same
                     config -> same dir -> idempotent, never clobbers a different
                     config. Reuses reproducibility-sketch/provenance.py.

  collate(...)    -> comparison. Walks runs/*/ and stacks every run into ONE
                     long-format table, each row tagged with that run's params,
                     so "compare backends / filters" is a group-by, not a manual
                     diff. `pivot_mean` reshapes to wide (chrom x run) for eyeballing.
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
from pathlib import Path

from config import CoverageConfig, DEFAULT_EXCLUDE
from result import ChromCoverage

# Reuse the existing provenance sidecar rather than reinventing it.
_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "dev" / "reproducibility-sketch"))
import provenance  # noqa: E402

# The config fields that DEFINE a run (i.e. change these -> a different output).
# Inputs are folded into the hash separately via their file identity.
RUN_PARAM_FIELDS = (
    "backend",
    "min_mapping_quality",
    "include_flags",
    "exclude_flags",
    "exclude_all_flags",
)

DATA_COLUMNS = ["chrom", "length", "bases", "mean_coverage"]


def run_params(config: CoverageConfig) -> dict:
    return {f: getattr(config, f) for f in RUN_PARAM_FIELDS}


def _stable_input_id(config: CoverageConfig) -> dict:
    # Size (not mtime) so the key is stable across a harmless `touch`; the full
    # provenance sidecar still records path/size/mtime/hash for real verification.
    return {
        "cram": {"name": config.cram.name, "size": config.cram.stat().st_size},
        "reference": {"name": config.reference.name, "size": config.reference.stat().st_size},
    }


def run_key(config: CoverageConfig, length: int = 8) -> str:
    """Short deterministic digest of what makes this run this run."""
    payload = {"params": run_params(config), "inputs": _stable_input_id(config)}
    blob = json.dumps(payload, sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()[:length]


def run_slug(config: CoverageConfig) -> str:
    """Human-scannable name: backend + only params that DIFFER from default,
    then the run_key hash as a suffix.

    The slug is friendly but lossy (it hides defaults and input identity); the
    hash suffix is what keeps the name unique + idempotent -- two runs that look
    identical in the slug but differ in a hidden param/input still get distinct
    dirs. Samtools-style letters (q/f/F/G) mirror the CLI flags; masks in hex.
    """
    parts = [config.backend]
    if config.min_mapping_quality:
        parts.append(f"q{config.min_mapping_quality}")
    if config.exclude_flags != DEFAULT_EXCLUDE:
        parts.append(f"F{config.exclude_flags:#x}")
    if config.include_flags:
        parts.append(f"f{config.include_flags:#x}")
    if config.exclude_all_flags:
        parts.append(f"G{config.exclude_all_flags:#x}")
    return "-".join(parts) + "-" + run_key(config)


def run_dirname(config: CoverageConfig, style: str = "slug") -> str:
    """Directory name for a run. 'slug' = readable + hash; 'hash' = bare digest."""
    if style == "hash":
        return run_key(config)
    if style == "slug":
        return run_slug(config)
    raise ValueError(f"unknown name style {style!r}; choose 'slug' or 'hash'")


def karyotypic_key(chrom: str):
    """Sort chr1..22, X, Y, M/MT, then unplaced/alt contigs alphabetically.

    Stable ordering makes even a raw `diff` of two per-run TSVs meaningful.
    """
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name.upper() in special:
        return (1, special[name.upper()], "")
    return (2, 0, chrom)  # alts / unplaced


def _header_lines(record: dict) -> list[str]:
    """VCF/samtools-style `#` comment header carrying run identity."""
    tool, code, params = record["tool"], record["code"], record["params"]
    kv = "  ".join(f"{k}={v}" for k, v in params.items())
    cram = record["inputs"]["cram"]
    ref = record["inputs"]["reference"]
    return [
        f"# {tool['name']} v{tool['version']}  commit={code['commit']}  dirty={code['dirty']}",
        f"# {kv}",
        f"# cram={cram['path']} ({cram['size_bytes']} bytes)",
        f"# reference={ref['path']}",
    ]


def write_run(
    rows: list[ChromCoverage],
    config: CoverageConfig,
    runs_dir: Path = Path("runs"),
    force: bool = False,
    name_style: str = "slug",
) -> Path:
    """Write runs/<name>/coverage.tsv + sidecar. Returns the run directory.

    `name_style` picks the directory name: 'slug' (readable + hash suffix) or
    'hash' (bare digest). Both end in the same hash, so if the directory already
    exists and force is False it's a no-op (same config already computed) -- the
    reproducibility payoff of a deterministic key.
    """
    out_dir = Path(runs_dir) / run_dirname(config, name_style)
    tsv = out_dir / "coverage.tsv"
    if tsv.exists() and not force:
        return out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    record = provenance.build_provenance(
        params=run_params(config),
        cram=config.cram,
        crai=config.index,
        reference=config.reference,
        outputs=[],  # filled after we know the file; sidecar re-stamped below
    )

    ordered = sorted(rows, key=lambda r: karyotypic_key(r.chrom))
    with tsv.open("w", newline="") as fh:
        for line in _header_lines(record):
            fh.write(line + "\n")
        w = csv.DictWriter(fh, fieldnames=DATA_COLUMNS, delimiter="\t")
        w.writeheader()
        for r in ordered:
            cells = r.as_row()
            w.writerow({c: cells[c] for c in DATA_COLUMNS})

    # Re-stamp provenance now that the output exists, and write the sidecar.
    record["outputs"] = [provenance.file_identity(tsv)]
    provenance.write_sidecar(record, tsv)
    return out_dir


# --- comparison ------------------------------------------------------------

def _read_run(run_dir: Path) -> list[dict]:
    """One run dir -> long-format rows tagged with its params (from the sidecar)."""
    tsv = run_dir / "coverage.tsv"
    sidecar = tsv.with_suffix(tsv.suffix + ".provenance.json")
    params = json.loads(sidecar.read_text())["params"] if sidecar.exists() else {}

    rows = []
    with tsv.open() as fh:
        reader = csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        for rec in reader:
            rows.append(
                {
                    "run_id": run_dir.name,
                    **params,
                    "chrom": rec["chrom"],
                    "length": int(rec["length"]),
                    "bases": int(rec["bases"]),
                    "mean_coverage": float(rec["mean_coverage"]),
                }
            )
    return rows


def collate(runs_dir: Path = Path("runs")) -> list[dict]:
    """Stack every run under runs_dir into one long-format table."""
    out: list[dict] = []
    for run_dir in sorted(Path(runs_dir).iterdir()):
        if (run_dir / "coverage.tsv").exists():
            out.extend(_read_run(run_dir))
    return out


def pivot_mean(long_rows: list[dict]) -> tuple[list[str], dict[str, dict[str, float]]]:
    """Reshape long rows to wide: chrom -> {run_id: mean_coverage}.

    Returns (run_ids, table) where table[chrom][run_id] = mean. Stdlib-only; a
    real tool would hand `long_rows` to pandas/polars for this.
    """
    run_ids = sorted({r["run_id"] for r in long_rows})
    table: dict[str, dict[str, float]] = {}
    for r in long_rows:
        table.setdefault(r["chrom"], {})[r["run_id"]] = r["mean_coverage"]
    ordered = {c: table[c] for c in sorted(table, key=karyotypic_key)}
    return run_ids, ordered
