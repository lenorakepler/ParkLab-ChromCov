"""
Data acquisition: download callability stratification files
Should also have ability to download CRAM files.
"""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

# Park Lab SMaHT_Regional_Categorization (GRCh38) callability tiers.
SMAHT_BASE = "https://raw.githubusercontent.com/parklab/SMaHT_Regional_Categorization/main"
SMAHT_STRATA = {
    "easy": "SMaHT_easy_hg38.bed.gz",
    "difficult": "SMaHT_difficult_hg38.bed.gz",
    "extreme": "SMaHT_extreme_hg38.bed.gz",
}

def download(url: str, dest: str | Path, force: bool = False) -> tuple[Path, bool]:
    """Download `url` to `dest` atomically. Returns (path, downloaded?) where
    downloaded is False if the file already existed and force is False."""
    dest = Path(dest)
    if dest.exists() and not force:
        return dest, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, tmp.open("wb") as fh:
        shutil.copyfileobj(resp, fh)
    tmp.replace(dest)
    return dest, True


def fetch_strata(dest_dir: str | Path = "data", force: bool = False) -> dict[str, Path]:
    """Download the three SMaHT strata BEDs into `dest_dir`. Returns
    {label: local path}. Prints one line per file (downloaded vs skipped)."""
    dest_dir = Path(dest_dir)
    out: dict[str, Path] = {}
    for label, name in SMAHT_STRATA.items():
        path, downloaded = download(f"{SMAHT_BASE}/{name}", dest_dir / name, force=force)
        size = path.stat().st_size
        status = "downloaded" if downloaded else "already present"
        print(f"[fetch] {label:9s} {status}: {path} ({size:,} bytes)")
        out[label] = path
    return out


def default_strata_paths(dest_dir: str | Path = "data") -> dict[str, Path]:
    """Expected local paths of the standard SMaHT strata BEDs (whether or not they
    exist yet). This is what `chromcov coverage --strata` resolves to, so the tiers
    never have to be named on the CLI -- it's the fixed easy/difficult/extreme set."""
    dest_dir = Path(dest_dir)
    return {label: dest_dir / name for label, name in SMAHT_STRATA.items()}
