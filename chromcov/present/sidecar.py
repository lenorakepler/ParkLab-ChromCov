"""
run.json provenance sidecar (extracted from cli.py helpers).

Writes a self-describing, re-runnable record for a --full run: tool + interpreter
version, git provenance, the resolved Config, contigs, baseline, and flagged
chromosomes -- so `chromcov plot` (and a human) can reconstruct the run from the
output dir alone.
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

RUN_SIDECAR = "run.json"


def _resolve(p: str | None) -> Path | None:
    return Path(p).expanduser().resolve() if p else None


def _tool_version() -> dict:
    """Package version + interpreter -- pins the analysis code in the sidecar."""
    try:
        pkg_version = version("parklab-chromcov")
    except PackageNotFoundError:
        pkg_version = None
    return {"name": "parklab-chromcov", "version": pkg_version, "python": sys.version.split()[0]}


def _run_git(args: list[str]) -> str | None:
    try:
        out = subprocess.run(["git", *args], cwd=Path(__file__).resolve().parent,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def _git_provenance() -> dict:
    """Exact commit + whether the tree was dirty. A commit SHA is meaningless if
    uncommitted changes were on disk, so `dirty` is the part that protects repro."""
    status = _run_git(["status", "--porcelain"])
    return {
        "commit": _run_git(["rev-parse", "HEAD"]),
        "branch": _run_git(["rev-parse", "--abbrev-ref", "HEAD"]),
        "describe": _run_git(["describe", "--tags", "--dirty", "--always"]),
        "dirty": bool(status) if status is not None else None,
    }


def write_run_sidecar(path: Path, result, config_file) -> None:
    """run.json for a --full run: self-describing + re-runnable. Embeds the
    resolved Config so `chromcov plot` (and a human) can reconstruct the run from
    the output dir alone. Takes a finalized RunResult."""
    value, source = result.baseline
    record = {
        "schema": "chromcov.run/3",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tool": _tool_version(),
        "code": _git_provenance(),
        "config_file": str(_resolve(config_file)) if config_file else None,
        "config": result.cfg.model_dump(mode="json"),
        "chromosomes": result.chroms,
        "baseline": {"value": round(value, 4), "source": source},
        "flagged": {c: fl for c, fl in result.flagged},
    }
    path.write_text(json.dumps(record, indent=2, sort_keys=True))
