"""
chromcov command-line interface (Click)
    
    # DEFAULT: per-chromosome mean table, no per-base output
    chromcov coverage --cram ... --reference ... -o ...
    (leave out -o for stdout)

    # EXTRAS: output per-base bedgraphs, stats, plots, stratification
              this is also resume-able given an existing directory
    chromcov coverage ... --full --outdir out/

    # RE-GRAPH: re-generate plots from existing bedgraphs
    chromcov plot --outdir out/

    # FETCH STRATIFICATION: download files that categorize the genome 
                            based on practical ability to variant call
                            (github.com/parklab/SMaHT_Regional_Categorization)
    chromcov fetch strata

Handles config loading (from Config.config), overriding, and outputting to sidecar
json for reproducibility. Dispatches to desired functions.
"""
from __future__ import annotations

import json
import subprocess
import sys
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import click
from pydantic import ValidationError

from . import fetch, report
from .config import Config
from .coverage import run_coverage
from .qc_report import QCReport
from .strata import Strata

# --------------------------------------------------- #
# ---------- config handling, CLI override ---------- #
# ----------- Read, Override, -> Sidecar ------------ #
# --------------------------------------------------- #

RUN_SIDECAR = "run.json"

def _resolve(p: str | None) -> Path | None:
    return Path(p).expanduser().resolve() if p else None

def _overrides(**kw) -> dict:
    """CLI options actually given -> a dict of Config field names. Paths resolved
    absolutely; None/empty dropped so model defaults (and any --config) apply."""
    ov: dict = {}
    for key in ("cram", "reference", "index"):
        if kw.get(key):
            ov[key] = str(_resolve(kw[key]))
    if kw.get("min_mapq") is not None:
        ov["min_mapping_quality"] = kw["min_mapq"]
    if kw.get("chroms"):
        ov["chroms"] = tuple(kw["chroms"].split(","))
    if kw.get("window"):
        ov["window"] = kw["window"]
    if kw.get("strata"):
        ov["strata"] = kw["strata"]     # already a {label: path} mapping (resolved in `coverage`)
    if kw.get("outdir"):
        ov["outdir"] = str(Path(kw["outdir"]).expanduser().resolve())
    return ov

def _load(config, **overrides) -> Config:
    try:
        return Config.load(config, overrides=_overrides(**overrides))
    except ValidationError as e:
        raise click.UsageError(
            "invalid or incomplete configuration -- provide --cram and --reference "
            f"(or a --config that sets them).\n{e}"
        )


def _emit(rows, output) -> None:
    """Write the coverage table: stdout by default (or when -o '-'), else to FILE."""
    frame = report.coverage_frame(rows)
    if output in (None, "-"):
        report.write_table(frame, "-")
        return
    target = Path(output).expanduser().resolve()
    report.write_table(frame, target)
    click.echo(f"wrote {target}", err=True)


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


def _write_run_sidecar(path: Path, report: QCReport, config_file) -> None:
    """run.json for a --full run: self-describing + re-runnable. Embeds the
    resolved Config so `chromcov plot` (and a human) can reconstruct the run from
    the output dir alone."""
    baseline, source = report.baseline()
    record = {
        "schema": "chromcov.run/3",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tool": _tool_version(),
        "code": _git_provenance(),
        "config_file": str(_resolve(config_file)) if config_file else None,
        "config": report.cfg.model_dump(mode="json"),
        "chromosomes": report.chroms,
        "baseline": {"value": round(baseline, 4), "source": source},
        "flagged": {c: fl for c, fl in report.flagged},
    }
    path.write_text(json.dumps(record, indent=2, sort_keys=True))


def input_options(func):
    func = click.option("--cram", type=click.Path(exists=True, dir_okay=False),
                        default=None, help="path to the CRAM")(func)
    func = click.option("--reference", type=click.Path(exists=True, dir_okay=False),
                        default=None, help="reference FASTA the CRAM was compressed against")(func)
    func = click.option("--index", type=click.Path(exists=True, dir_okay=False),
                        default=None, help="CRAM index (default: <cram>.crai)")(func)
    func = click.option("--min-mapq", "min_mapq", type=int, default=None,
                        help="minimum mapping quality (-Q)")(func)
    func = click.option("--config", type=click.Path(exists=True, dir_okay=False),
                        default=None, help="run-config YAML (CLI options override its values)")(func)
    return func

# --------------------------------------------------- #
# ------------------ Dispatching -------------------- #
# --------------------------------------------------- #

@click.group()
@click.version_option(package_name="parklab-chromcov", prog_name="chromcov")
def main() -> None:
    """Per-chromosome average coverage from a CRAM, with optional QC extensions."""


@main.command()
@input_options
@click.option("--chroms", default=None, help="comma-separated contig subset (default: all)")
@click.option("--full", is_flag=True,
              help="beyond the mean: per-base depth once (as bedgraphs, which also "
                   "checkpoint resume) -> stats + windows + strata + copy number + plots")
@click.option("--window", type=int, default=None, help="[--full] windowed-mean bin size (bp)")
@click.option("--strata", is_flag=True,
              help="[--full] stratify by the SMaHT easy/difficult/extreme callability "
                   "tiers: colors the scatter by tier and uses the 'easy' tier as the "
                   "copy-number baseline. Uses the BEDs in --strata-dir (get them with "
                   "`chromcov fetch strata`) -- no tier names needed.")
@click.option("--strata-dir", "strata_dir", type=click.Path(file_okay=False), default="data",
              help="directory holding the SMaHT strata BEDs for --strata (default: ./data)")
@click.option("--jobs", "-j", type=int, default=1,
              help="compute contigs in parallel across N processes (default 1)")
@click.option("--force", is_flag=True, help="[--full] recompute bedgraphs even if present")
@click.option("--outdir", type=click.Path(file_okay=False), default=None,
              help="[--full] output root (default: ./out); bedgraphs -> <outdir>/perbase/")
@click.option("--output", "-o", default=None,
              help="write the mean table to FILE ('-' for stdout; the default)")
def coverage(cram, reference, index, min_mapq, config, chroms, full, window, strata,
             strata_dir, jobs, force, outdir, output) -> None:
    """Per-chromosome average coverage.

    By default reports the mean per chromosome (bases / length) to stdout. `--full`
    additionally computes per-base depth and derives robust stats, callability
    strata, approximate copy number, QC flags, and plots -- writing per-base
    bedgraphs under <outdir>/perbase/ that make the run resumable and the graphs
    incremental. The mean is identical to the default run."""
    # --strata is all-or-nothing: on -> the fixed SMaHT tier set (resolved from
    # --strata-dir); off -> no stratification (unless a --config supplies its own).
    strata_map = None
    if strata:
        paths = fetch.default_strata_paths(strata_dir)
        missing = [p.name for p in paths.values() if not p.exists()]
        if missing:
            raise click.UsageError(
                f"--strata needs the SMaHT callability BEDs in {strata_dir}/ (missing: "
                f"{', '.join(missing)}). Run `chromcov fetch strata` (or set --strata-dir).")
        strata_map = {label: str(p) for label, p in paths.items()}

    cfg = _load(config, cram=cram, reference=reference, index=index, min_mapq=min_mapq,
                chroms=chroms, window=window, strata=strata_map, outdir=outdir)

    if not full:
        _emit(run_coverage(cfg, jobs=jobs), output)
        return

    report = QCReport(cfg, Strata.from_arg(cfg.strata))
    bedgraph_dir = Path(cfg.outdir) / "perbase"
    pf = report.run(bedgraph_dir=bedgraph_dir, jobs=jobs, force=force)
    click.echo(f"[preflight] {pf['reference_check']['status']}", err=True)

    outdir_p = Path(cfg.outdir)
    written = report.write_outputs(outdir_p)
    _write_run_sidecar(outdir_p / RUN_SIDECAR, report, config)
    if output:
        _emit(report.coverage_rows(), output)

    for line in report.summary_lines():
        click.echo(line, err=True)
    click.echo(f"per-base bedgraphs: {bedgraph_dir}", err=True)
    click.echo(f"outputs: {outdir_p}", err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)


@main.command()
@click.option("--outdir", type=click.Path(file_okay=False), default="out",
              help="a --full run dir (holds run.json + perbase/); default ./out")
def plot(outdir) -> None:
    """(Re)build the tables and plots from the bedgraphs already under --outdir.

    Reads the run's config from run.json and reduces whatever contigs have a
    bedgraph in perbase/ -- so running more contigs (coverage --full ...) and then
    `plot` updates the graphs, with the copy-number baseline recomputed over every
    contig present. No CRAM needed."""
    outdir_p = Path(outdir)
    sidecar = outdir_p / RUN_SIDECAR
    if not sidecar.exists():
        raise click.UsageError(
            f"no {RUN_SIDECAR} under {outdir_p}/ -- run `chromcov coverage --full "
            f"--outdir {outdir}` first")
    cfg = Config.model_validate(json.loads(sidecar.read_text())["config"])

    report = QCReport(cfg, Strata.from_arg(cfg.strata))
    report.gather(outdir_p / "perbase")
    if not report.chroms:
        raise click.UsageError(f"no per-base bedgraphs under {outdir_p / 'perbase'}/")
    written = report.write_outputs(outdir_p)
    for line in report.summary_lines():
        click.echo(line, err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)


@main.group(name="fetch")
def fetch_group() -> None:
    """Download inputs a clean clone needs (callability strata, ...)."""


@fetch_group.command(name="strata")
@click.option("--dest", type=click.Path(file_okay=False), default="data",
              help="directory to download into (default: ./data)")
@click.option("--force", is_flag=True, help="re-download even if the files exist")
def fetch_strata_cmd(dest, force) -> None:
    """Download the Park Lab SMaHT easy/difficult/extreme hg38 BEDs."""
    fetch.fetch_strata(dest, force=force)
    hint = "  chromcov coverage ... --full --strata"
    if str(dest) not in ("data", "./data"):
        hint += f" --strata-dir {dest}"
    click.echo("\nuse them with:\n" + hint)


if __name__ == "__main__":
    main()
