"""
Coverage backends: one shared `CoverageConfig` in, one normalized
`list[ChromCoverage]` out. Two implementations behind a common ABC so downstream
code is backend-agnostic and "run both, diff the means" is a real cross-check:

  NativeBackend    thin adapter over the hand-rolled pysam calculator
                   (chromcov.read_filter.calc_cov).
  MosdepthBackend  subprocess wrapper over the `mosdepth` binary on PATH + its
                   .mosdepth.summary.txt parser.

The `run(config) -> list[ChromCoverage]` contract is the abstract method. Where a
config knob can't be honored, the backend RAISES rather than silently diverging
(silent divergence defeats the whole point of having two backends).
"""
from __future__ import annotations

import csv
import shutil
import subprocess
import tempfile
from abc import ABC, abstractmethod
from pathlib import Path

import pysam

from .config import CoverageConfig
from .read_filter import ReadFilter, calc_cov
from .result import ChromCoverage


class CoverageBackend(ABC):
    """A coverage engine: config in, normalized per-chromosome rows out."""

    name: str

    @abstractmethod
    def run(self, config: CoverageConfig) -> list[ChromCoverage]:
        ...


class NativeBackend(CoverageBackend):
    """Hand-rolled pysam calculator. Translates the shared config into a
    ReadFilter once, opens the CRAM once, loops chromosomes. All the actual
    coverage math stays in `calc_cov`; this only bridges config <-> that API."""

    name = "native"

    def run(self, config: CoverageConfig) -> list[ChromCoverage]:
        rf = ReadFilter(
            include_flags=config.include_flags,
            exclude_flags=config.exclude_flags,
            exclude_all_flags=config.exclude_all_flags,
            min_mapping_quality=config.min_mapping_quality,
        )

        cram = pysam.AlignmentFile(
            str(config.cram), "rc",
            reference_filename=str(config.reference),
            index_filename=str(config.index),
        )

        chroms = config.chroms if config.chroms is not None else cram.references

        rows: list[ChromCoverage] = []
        for chrom in chroms:
            base_depth, total_depth, _ = calc_cov(cram, chrom, rf, per_base=config.per_base)
            rows.append(
                ChromCoverage(
                    chrom=chrom,
                    length=cram.get_reference_length(chrom),
                    bases=int(total_depth),
                    backend=self.name,
                    base_depth=base_depth,
                )
            )
        return rows


class MosdepthBackend(CoverageBackend):
    """Shell out to the `mosdepth` binary on PATH (never `docker run`) and parse
    its summary. Reproducibility of the binary belongs to the environment layer
    (the pinned Docker/conda env the CWL runner provides), keeping this symmetric
    with the native backend (both in-process, same env)."""

    name = "mosdepth"

    def run(self, config: CoverageConfig) -> list[ChromCoverage]:
        if shutil.which("mosdepth") is None:
            raise RuntimeError(
                "mosdepth not found on PATH. Install it in the environment "
                "(conda: `mamba install -c bioconda mosdepth`, or the biocontainer/Docker image)."
            )

        # mosdepth writes several files next to <prefix>; isolate them in a temp
        # dir so the caller only deals with the parsed rows.
        with tempfile.TemporaryDirectory() as td:
            prefix = Path(td) / "cov"
            subprocess.run(self._build_argv(config, prefix), check=True)
            return self._parse_summary(prefix.with_suffix(".mosdepth.summary.txt"))

    def _build_argv(self, config: CoverageConfig, prefix: Path) -> list[str]:
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
            # mosdepth -c restricts to a single contig. A multi-contig subset needs
            # a BED via --by; keep it honest by raising rather than silently
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

    def _parse_summary(self, summary_path: Path) -> list[ChromCoverage]:
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
                        backend=self.name,
                    )
                )
        return rows


# Registry: name -> singleton backend instance.
BACKENDS: dict[str, CoverageBackend] = {
    b.name: b for b in (NativeBackend(), MosdepthBackend())
}


def get_backend(name: str) -> CoverageBackend:
    try:
        return BACKENDS[name]
    except KeyError:
        raise ValueError(f"unknown backend {name!r}; choose from {sorted(BACKENDS)}")
