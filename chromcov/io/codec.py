"""
Reference checksum / verification (from validate.py).

Confirms the supplied reference is the one the CRAM was compressed against, by
comparing the CRAM's @SQ M5 tags to the reference's per-sequence MD5 (from a
`.dict` sidecar, or computed from the .fa). A mismatch means coverage would be
silently wrong, so it raises `PreflightError`.
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import pysam

from . import alignment


class PreflightError(Exception):
    """A hard input problem that must stop the run."""


def _reference_m5_from_dict(reference: Path) -> dict[str, str]:
    """Per-sequence M5 from a samtools .dict sidecar, if present (no hashing)."""
    for cand in (reference.with_suffix(".dict"), Path(str(reference) + ".dict")):
        if cand.exists():
            out = {}
            for line in cand.read_text().splitlines():
                if line.startswith("@SQ"):
                    f = dict(kv.split(":", 1) for kv in line.split("\t")[1:] if ":" in kv)
                    if "SN" in f and "M5" in f:
                        out[f["SN"]] = f["M5"]
            return out

    return {}


def _reference_m5_computed(reference: Path, names: list[str]) -> dict[str, str]:
    """Compute M5 for the named sequences by reading the .fa (one pass)."""
    fa = pysam.FastaFile(str(reference))
    out = {}
    for name in names:
        if name in fa.references:
            seq = fa.fetch(name).upper().encode()
            out[name] = hashlib.md5(seq).hexdigest()
    fa.close()
    return out


def check_reference(config, verify: str = "auto") -> dict:
    """Compare CRAM @SQ M5 to the reference. verify: 'auto' (dict only, warn if
    absent), 'full' (hash the .fa if no dict), 'skip'."""
    sq = alignment.read_header(config).get("SQ", [])
    cram_m5 = {s["SN"]: s.get("M5") for s in sq}
    if verify == "skip":
        return {"status": "skipped", "checked": 0}

    ref_m5 = _reference_m5_from_dict(config.reference)
    source = "dict"
    if not ref_m5:
        if verify == "full":
            ref_m5 = _reference_m5_computed(config.reference, list(cram_m5))
            source = "computed"
        else:
            return {"status": "warn",
                    "message": "no reference .dict found; run `samtools dict` or pass "
                               "verify='full' to hash the .fa for M5 verification"}

    mismatches = [sn for sn, m5 in cram_m5.items()
                  if m5 and sn in ref_m5 and ref_m5[sn] != m5]
    if mismatches:
        raise PreflightError(
            f"reference does not match the CRAM: M5 mismatch on {len(mismatches)} "
            f"sequence(s), e.g. {mismatches[:3]}. The supplied reference is not the "
            "one this CRAM was compressed against; coverage would be silently wrong."
        )
    missing = [sn for sn, m5 in cram_m5.items() if not m5]
    return {"status": "ok", "source": source,
            "checked": sum(1 for sn, m5 in cram_m5.items() if m5 and sn in ref_m5),
            "missing_m5": missing}
