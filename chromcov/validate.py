"""
Input preflight validation:

Checks:
  1. inputs exist                (cram, reference, index)
  2. CRAM is coordinate-sorted   (@HD SO:coordinate)
  3. index is present            (.crai -- required by fetch(chrom))
  4. reference matches the CRAM  (@SQ M5 tags vs the reference's per-sequence MD5)
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pysam


class PreflightError(Exception):
    """A hard input problem that must stop the run."""

def _cram_header(config) -> dict:
    af = pysam.AlignmentFile(str(config.cram), "rc",
                             reference_filename=str(config.reference),
                             index_filename=str(config.index))
    hdr = af.header.to_dict()
    af.close()
    
    return hdr

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
    sq = _cram_header(config).get("SQ", [])
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

    hdr = _cram_header(config)
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
