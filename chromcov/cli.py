"""
chromcov command-line interface (Click)

    # DEFAULT: per-chromosome mean table, no per-base output
    # ---------------------------------------------------------------
    chromcov coverage --cram ... --reference ... -o ...
    (leave out -o for stdout)

    # EXTRAS: output per-base tracks, stats, plots, stratification
              this is also resume-able given an existing directory
    # ---------------------------------------------------------------
    chromcov coverage ... --full --outdir out/

    # RE-GRAPH: re-generate plots from existing tracks
    # ---------------------------------------------------------------
    chromcov plot --outdir out/

    # FETCH STRATIFICATION: download files that categorize the genome
    # ---------------------------------------------------------------
                            based on practical ability to variant call
                            (github.com/parklab/SMaHT_Regional_Categorization)
    chromcov fetch strata

Thin dispatch: parse -> Config.load -> pipeline.run -> present. Config loading,
CLI override, and the provenance sidecar are delegated to their modules.
"""
from __future__ import annotations

import json
from pathlib import Path

import click
from pydantic import ValidationError

from .categories import Strata
from .config.schema import Config
from .io import fetch
from .pipeline import Depth, Source, run as pipeline_run
from .present import frames, sidecar
from .present.sidecar import RUN_SIDECAR, _resolve

# ==============================================================================
# CONFIG HANDLING
# ==============================================================================
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

# ==============================================================================
# SETUP FUNCS
# ==============================================================================
def _emit(rows, output) -> None:
    """Write the coverage table: stdout by default (or when -o '-'), else to FILE."""
    frame = frames.coverage_frame(rows)
    if output in (None, "-"):
        frames.write_table(frame, "-")
        return
    target = Path(output).expanduser().resolve()
    frames.write_table(frame, target)
    click.echo(f"wrote {target}", err=True)


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

# ==============================================================================
# COVERAGE CALC DISPATCHING
# ==============================================================================
@click.group()
@click.version_option(package_name="parklab-chromcov", prog_name="chromcov")
def main() -> None:
    """Per-chromosome average coverage from a CRAM, with optional QC extensions."""

@main.command()
@input_options
@click.option("--chroms", default=None, help="comma-separated contig subset (default: all)")
@click.option("--full", is_flag=True,
              help="beyond the mean: per-base depth once (as tracks, which also "
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
@click.option("--force", is_flag=True, help="[--full] recompute tracks even if present")
@click.option("--outdir", type=click.Path(file_okay=False), default=None,
              help="[--full] output root (default: ./out); tracks -> <outdir>/perbase/")
@click.option("--output", "-o", default=None,
              help="write the mean table to FILE ('-' for stdout; the default)")
def coverage(cram, reference, index, min_mapq, config, chroms, full, window, strata,
             strata_dir, jobs, force, outdir, output) -> None:
    """
    Per-chromosome average coverage.

    By default reports the mean per chromosome (bases / length) to stdout.

    `--full` also computes per-base depth and derives robust stats, callability
    strata, approximate copy number, QC flags, and plots -- writing per-base tracks
    under <outdir>/perbase/ that make the run resumable and the graphs incremental.
    The mean is identical to the default run.
    """

    # ------------------------------------------------------------------------------
    # Stratify
    # ------------------------------------------------------------------------------
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

    # ------------------------------------------------------------------------------
    # Mean-only return
    # ------------------------------------------------------------------------------
    if not full:
        result = pipeline_run(cfg, depth=Depth.MEAN, jobs=jobs)
        _emit(result.coverage_rows(), output)
        return

    # ------------------------------------------------------------------------------
    # Full
    # ------------------------------------------------------------------------------
    bedgraph_dir = Path(cfg.outdir) / "perbase"
    result = pipeline_run(
        cfg, depth=Depth.FULL, 
        source=Source.ALIGNMENT,
        jobs=jobs,
        force=force, 
        bedgraph_dir=bedgraph_dir, 
        categories=Strata.from_arg(cfg.strata)
        )
    click.echo(f"[preflight] {result.preflight['reference_check']['status']}", err=True)

    # Output
    outdir_p = Path(cfg.outdir)
    written = frames.write_outputs(result, outdir_p)
    sidecar.write_run_sidecar(outdir_p / RUN_SIDECAR, result, config)
    if output:
        _emit(result.coverage_rows(), output)

    for line in frames.summary_lines(result):
        click.echo(line, err=True)

    click.echo(f"per-base tracks: {bedgraph_dir}", err=True)
    click.echo(f"outputs: {outdir_p}", err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)

# ==============================================================================
# PLOTTING
# ==============================================================================
@main.command()
@click.option("--outdir", type=click.Path(file_okay=False), default="out",
              help="a --full run dir (holds run.json + perbase/); default ./out")
def plot(outdir) -> None:
    """
    (Re)build the tables and plots from the tracks already under --outdir.

    Reads the run's config from run.json and reduces whatever contigs have a track
    in perbase/ -- so running more contigs (coverage --full ...) and then `plot`
    updates the graphs, with the copy-number baseline recomputed over every contig
    present. No CRAM needed.
    """
    outdir_p = Path(outdir)
    sidecar_path = outdir_p / RUN_SIDECAR
    if not sidecar_path.exists():
        raise click.UsageError(
            f"no {RUN_SIDECAR} under {outdir_p}/ -- run `chromcov coverage --full "
            f"--outdir {outdir}` first")
    cfg = Config.model_validate(json.loads(sidecar_path.read_text())["config"])

    result = pipeline_run(cfg, depth=Depth.FULL, source=Source.TRACKS,
                          bedgraph_dir=outdir_p / "perbase", categories=Strata.from_arg(cfg.strata))
    if not result.chroms:
        raise click.UsageError(f"no per-base tracks under {outdir_p / 'perbase'}/")
    written = frames.write_outputs(result, outdir_p)
    for line in frames.summary_lines(result):
        click.echo(line, err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)

# ==============================================================================
# CONFIG GENERATION
# ==============================================================================
@main.command(name="gen-config")
@click.option("--output", "-o", "output", type=click.Path(dir_okay=False),
              default="config.yaml", help="path to write (default: ./config.yaml)")
@click.option("--force", is_flag=True, help="overwrite if the file already exists")
def gen_config_cmd(output, force) -> None:
    """
    Write an editable run-config YAML pre-filled with every option at its
    current default.

    Values are read live from the Config model, so the generated file always
    matches the code's defaults -- edit it and pass it back with `--config`.
    """
    from .config.template import write_default_config

    path = Path(output)
    if path.exists() and not force:
        raise click.UsageError(f"{path} already exists; pass --force to overwrite")
    write_default_config(path)
    click.echo(f"wrote {path}")
    click.echo(f"edit it, then: chromcov coverage --config {path}", err=True)

# ==============================================================================
# STRATA FETCHING
# ==============================================================================
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
