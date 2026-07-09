"""
Input preflight validation (was validate.py, minus the checksum).

Checks:
  1. inputs exist                (cram, reference, index)
  2. CRAM is coordinate-sorted   (@HD SO:coordinate)
  3. index is present            (.crai -- required by fetch(chrom))
  4. reference matches the CRAM  -> delegated to io.codec.check_reference
"""
from __future__ import annotations

import sys
from pathlib import Path

from .io import alignment
from .io.codec import PreflightError, check_reference

__all__ = ["PreflightError", "preflight", "check_reference"]


def preflight(config, verify_reference: str | None = None) -> dict:
    """Run all checks. Raises PreflightError on hard failures; returns a report
    (suitable for the provenance sidecar) on success/soft-warnings.

    `verify_reference` (auto|full|skip) overrides the mode; None uses the run's
    `config.verify_reference`."""
    verify = verify_reference if verify_reference is not None else getattr(config, "verify_reference", "auto")
    for label, p in (("cram", config.cram), ("reference", config.reference), ("index", config.index)):
        if not Path(p).exists():
            hint = " (build with `samtools index`)" if label == "index" else ""
            raise PreflightError(f"{label} not found: {p}{hint}")

    hdr = alignment.read_header(config)
    so = hdr.get("HD", {}).get("SO")
    if so != "coordinate":
        raise PreflightError(
            f"CRAM is not coordinate-sorted (@HD SO={so!r}); fetch-by-region and "
            "the .crai require it. Run `samtools sort` (then re-index)."
        )

    reference = check_reference(config, verify)
    report = {"sorted": True, "index": str(config.index), "reference_check": reference}
    if reference["status"] == "warn":
        # stderr so the coverage table on stdout stays clean/pipeable.
        print(f"[preflight] WARNING: {reference['message']}", file=sys.stderr)

    return report
