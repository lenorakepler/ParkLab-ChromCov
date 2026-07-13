"""
Data acquisition: download the Park Lab inputs a clean clone needs -- the
COLO829T test CRAM (+ index) and its GRCh38 reference, plus the callability
strata BEDs. Everything lands in `data/` by default, which is where the model
defaults look for it.
"""
from __future__ import annotations

import shutil
import urllib.request
from pathlib import Path

# Where downloads land, and where the Config model's input defaults point.
DATA_DIR = "data"

# Park Lab COLO829T interview inputs (GRCh38), hosted on S3. These are the files
# `Config` defaults to (data/<name>), so a fetch-then-run needs no flags.
PARKLAB_BASE = "https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW"
PARKLAB_INPUTS = {
    "cram": "COLO829T_TEST.cram",
    "index": "COLO829T_TEST.cram.crai",
    "reference": "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
}

# Park Lab SMaHT_Regional_Categorization (GRCh38) callability tiers.
SMAHT_BASE = "https://raw.githubusercontent.com/parklab/SMaHT_Regional_Categorization/main"
SMAHT_STRATA = {
    "easy": "SMaHT_easy_hg38.bed.gz",
    "difficult": "SMaHT_difficult_hg38.bed.gz",
    "extreme": "SMaHT_extreme_hg38.bed.gz",
}

def download(url: str, dest: str | Path, force: bool = False) -> tuple[Path, bool]:
    """
    Download `url` to `dest` atomically. Returns (path, downloaded?) where
    downloaded is False if the file already existed and force is False.
    """
    dest = Path(dest)
    if dest.exists() and not force:
        return dest, False
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url) as resp, tmp.open("wb") as fh:
        shutil.copyfileobj(resp, fh)
    tmp.replace(dest)
    return dest, True


def default_input_paths(dest_dir: str | Path = DATA_DIR) -> dict[str, Path]:
    """
    Expected local paths of the Park Lab COLO829T inputs (whether or not they
    exist yet). These are exactly the paths the Config model defaults to, so
    `chromcov fetch inputs` followed by a bare `chromcov coverage` just works.
    """
    dest_dir = Path(dest_dir)
    return {label: dest_dir / name for label, name in PARKLAB_INPUTS.items()}


def fetch_inputs(dest_dir: str | Path = DATA_DIR, force: bool = False) -> dict[str, Path]:
    """
    Download the Park Lab COLO829T test CRAM, its .crai index, and the GRCh38
    reference FASTA into `dest_dir`. Returns {label: local path}. Prints one line
    per file (downloaded vs skipped). NOTE: the CRAM + reference are large (tens of
    GB); expect this to take a while.
    """
    dest_dir = Path(dest_dir)
    out: dict[str, Path] = {}
    for label, name in PARKLAB_INPUTS.items():
        path, downloaded = download(f"{PARKLAB_BASE}/{name}", dest_dir / name, force=force)
        size = path.stat().st_size
        status = "downloaded" if downloaded else "already present"
        print(f"[fetch] {label:9s} {status}: {path} ({size:,} bytes)")
        out[label] = path
    return out

def fetch_strata(dest_dir: str | Path = "data", force: bool = False) -> dict[str, Path]:
    """
    Download the three SMaHT strata BEDs into `dest_dir`. Returns
    {label: local path}. Prints one line per file (downloaded vs skipped).
    """
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
    """
    Expected local paths of the standard SMaHT strata BEDs (whether or not they
    exist yet). This is what `chromcov coverage --strata` resolves to, so the tiers
    never have to be named on the CLI -- it's the fixed easy/difficult/extreme set.
    """
    dest_dir = Path(dest_dir)
    return {label: dest_dir / name for label, name in SMAHT_STRATA.items()}
