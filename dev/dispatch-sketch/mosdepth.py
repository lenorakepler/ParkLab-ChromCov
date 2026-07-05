"""
mosdepth backend: shell out to the `mosdepth` binary and parse its summary.

Design choice: we call `mosdepth` as a subprocess against a binary on PATH; we
do NOT shell out to `docker run`. Reproducibility of the binary belongs to the
environment layer (the pinned Docker/conda env that the CWL runner provides),
not to application code. That keeps this backend symmetric with the native one
(both run in-process, same env) and avoids docker-in-docker once the whole tool
is itself containerized by the workflow engine.

We use `--no-per-base` unless the config asks for per-base output: for the
"average per chromosome" answer, the summary file is all we need and is far
faster (mosdepth doesn't write the big per-base BED).
"""
from __future__ import annotations
import csv
import shutil
import subprocess
import tempfile
from pathlib import Path

from config import CoverageConfig
from result import ChromCoverage


def _build_argv(config: CoverageConfig, prefix: Path) -> list[str]:
    if config.exclude_all_flags:
        # mosdepth has no -G / "exclude only if ALL bits set" equivalent. Fail
        # loudly rather than quietly producing numbers that differ from native.
        raise ValueError(
            "mosdepth backend cannot honor exclude_all_flags (-G); "
            "unset it or use backend='native'."
        )

    argv = ["mosdepth"]
    if not config.per_base:
        argv.append("--no-per-base")
    argv += [
        "--threads", str(config.threads),
        "--fasta", str(config.reference),   # CRAM decode needs the reference
        "--mapq", str(config.min_mapping_quality),
        "--flag", str(config.exclude_flags),
    ]
    if config.chroms is not None:
        # mosdepth -c restricts to a single contig. A multi-contig subset needs a
        # BED via --by; keep the sketch honest by raising rather than silently
        # computing the whole genome.
        if len(config.chroms) != 1:
            raise NotImplementedError(
                "mosdepth backend supports a single-contig subset (-c) only; "
                f"got {config.chroms}. Use a BED (--by) or run all contigs."
            )
        argv += ["--chrom", config.chroms[0]]
    if config.include_flags:
        argv += ["--include-flag", str(config.include_flags)]
    argv += [str(prefix), str(config.cram)]
    return argv


def _parse_summary(summary_path: Path, backend: str) -> list[ChromCoverage]:
    """Parse <prefix>.mosdepth.summary.txt.

    Columns: chrom, length, bases, mean, min, max. Skip the trailing `total`
    aggregate row and any `_region` rows (only present when a BED is passed).
    We recompute `mean` from bases/length in ChromCoverage rather than trusting
    mosdepth's rounded column, so it matches the native backend exactly.
    """
    rows: list[ChromCoverage] = []
    with open(summary_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for rec in reader:
            chrom = rec["chrom"]
            if chrom == "total" or chrom.endswith("_region"):
                continue
            rows.append(
                ChromCoverage(
                    chrom=chrom,
                    length=int(rec["length"]),
                    bases=int(rec["bases"]),
                    backend=backend,
                )
            )
    return rows


def run(config: CoverageConfig) -> list[ChromCoverage]:
    if shutil.which("mosdepth") is None:
        raise RuntimeError(
            "mosdepth not found on PATH. Install it in the environment "
            "(conda: `mamba install -c bioconda mosdepth`, or the biocontainer/Docker image)."
        )

    # mosdepth writes several files next to <prefix>; isolate them in a temp dir
    # so the caller only deals with the parsed rows.
    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / "cov"
        subprocess.run(_build_argv(config, prefix), check=True)
        return _parse_summary(prefix.with_suffix(".mosdepth.summary.txt"), config.backend)
