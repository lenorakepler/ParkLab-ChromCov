"""
provenance.py — capture "what code + what params + what inputs" alongside every
coverage result, so any output TSV can be traced back to an exact, reproducible run.

Design notes (calibrated to this assignment):
  * The CRAM is ~17 GB. Hashing it on every run is wasteful, so we record cheap
    *identity* (path + size + mtime) for big inputs and only sha256 the small
    files (the .crai, the .fai). Set hash_big=True to force full hashing.
  * CRAM is reference-based: the file cannot be decoded — and coverage cannot be
    reproduced — without the exact reference it was compressed against. Each CRAM
    @SQ line carries an M5 (MD5 of the uppercased, base-only sequence) tag. We
    extract those and (if pysam is available) let the caller verify the supplied
    reference matches, which is the single most important reproducibility check
    for CRAM-based coverage.

Stdlib-only by default; pysam is imported lazily and is optional.
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# Files below this size are cheap enough to hash on every run.
SMALL_FILE_BYTES = 256 * 1024 * 1024  # 256 MB


def _run_git(args: list[str], repo: Path) -> str | None:
    try:
        out = subprocess.run(
            ["git", *args],
            cwd=repo,
            capture_output=True,
            text=True,
            check=True,
        )
        return out.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None


def git_provenance(repo: Path | None = None) -> dict:
    """Exact commit + whether the working tree was dirty at run time.

    A commit SHA is meaningless if uncommitted changes were on disk, so the
    `dirty` flag is the part that actually protects reproducibility.
    """
    repo = repo or Path(__file__).resolve().parent
    commit = _run_git(["rev-parse", "HEAD"], repo)
    status = _run_git(["status", "--porcelain"], repo)
    branch = _run_git(["rev-parse", "--abbrev-ref", "HEAD"], repo)
    # Prefer a tag-derived version if one exists (mirrors `dunamai`/`git describe`).
    describe = _run_git(["describe", "--tags", "--dirty", "--always"], repo)
    return {
        "commit": commit,
        "branch": branch,
        "describe": describe,
        "dirty": bool(status) if status is not None else None,
    }


def tool_version() -> dict:
    """Package version + interpreter, so the analysis code itself is pinned."""
    version = None
    try:
        from importlib.metadata import version as _v

        version = _v("parklab-chromcov")
    except Exception:
        pass
    return {
        "name": "parklab-chromcov",
        "version": version,
        "python": sys.version.split()[0],
    }


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def file_identity(path: str | Path, hash_big: bool = False) -> dict:
    """Cheap, always-safe identity for an input; full sha256 only when small."""
    p = Path(path)
    st = p.stat()
    rec = {
        "path": str(p),
        "size_bytes": st.st_size,
        "mtime_utc": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
        "sha256": None,
    }
    if hash_big or st.st_size <= SMALL_FILE_BYTES:
        rec["sha256"] = _sha256(p)
    else:
        rec["sha256_skipped"] = f"file > {SMALL_FILE_BYTES} bytes; use hash_big=True"
    return rec


def cram_reference_ids(cram_path: str | Path, reference_path: str | Path) -> dict:
    """Per-sequence M5 tags from the CRAM header — the reference fingerprint.

    Returns the M5 (+UR/SP where present) for each @SQ line. Comparing these to
    the MD5s of the supplied reference sequences proves the reference matches the
    one the CRAM was written against. pysam is optional; absence is recorded, not
    fatal, so provenance capture never blocks a run.
    """
    try:
        import pysam
    except ImportError:
        return {"available": False, "reason": "pysam not importable"}

    af = pysam.AlignmentFile(
        str(cram_path), "rc", reference_filename=str(reference_path)
    )
    seqs = [
        {"name": sq["SN"], "length": sq["LN"], "M5": sq.get("M5"), "UR": sq.get("UR")}
        for sq in af.header.to_dict().get("SQ", [])
    ]
    af.close()
    return {"available": True, "reference_path": str(reference_path), "sequences": seqs}


def build_provenance(
    *,
    params: dict,
    cram: str | Path,
    crai: str | Path,
    reference: str | Path,
    outputs: list[str | Path],
    hash_big: bool = False,
) -> dict:
    """Assemble the full provenance record for one coverage run."""
    return {
        "schema": "chromcov.provenance/1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tool": tool_version(),
        "code": git_provenance(),
        "params": params,
        "inputs": {
            "cram": file_identity(cram, hash_big=hash_big),
            "crai": file_identity(crai, hash_big=hash_big),
            "reference": file_identity(reference, hash_big=hash_big),
        },
        "reference_verification": cram_reference_ids(cram, reference),
        "outputs": [file_identity(o, hash_big=hash_big) for o in outputs],
    }


def write_sidecar(record: dict, output_path: str | Path) -> Path:
    """Write `<output>.provenance.json` next to the result it describes."""
    p = Path(output_path)
    sidecar = p.with_suffix(p.suffix + ".provenance.json")
    sidecar.write_text(json.dumps(record, indent=2, sort_keys=True))
    return sidecar


if __name__ == "__main__":
    # Smoke test: emit provenance for this repo's own state (no CRAM needed).
    print(json.dumps({"tool": tool_version(), "code": git_provenance()}, indent=2))
