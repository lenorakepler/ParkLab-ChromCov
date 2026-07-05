"""
chromcov command-line interface (Click).

    chromcov coverage  --cram ... --reference ...      # fast per-chromosome mean (the deliverable)
    chromcov analyze   --cram ... --reference ...      # full QC suite (superset of coverage)
    chromcov collate                                   # compare archived runs

Path handling follows dev/output-location-conventions.md: inputs are resolved
absolutely (expanduser().resolve()); the coverage table defaults to stdout so it
is pipeable, or `--output FILE`; analysis outputs go under `--outdir` (default
./out). Nothing is written to the install directory.

`coverage` is what the CWL tool (dev/reproducibility-sketch/coverage.cwl) calls
as `baseCommand: [chromcov, coverage]` with `--cram --reference --min-mapq
--output`, so those options are kept stable.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import click
from pydantic import ValidationError

from . import dispatch, fetch, provenance
from .config import RunConfig
from .output import RunStore
from .perbase import PerBaseStore, analysis_key, analysis_slug, build_track
from .pipeline import CoverageAnalysis
from .result import TSV_COLUMNS, write_tsv
from .strata import Strata


def _resolve(p: str | None) -> Path | None:
    return Path(p).expanduser().resolve() if p else None


# The CLI only ever contributes *overrides*: dicts of the options actually given,
# mapped to config-field names. Defaults + validation live in the models; the
# config file (if any) is the base. This is the whole "config authoritative, CLI
# overrides" rule, and it lives only here.

def _coverage_overrides(*, cram=None, reference=None, index=None, min_mapq=None,
                        chroms=None, per_base=False) -> dict:
    ov: dict = {}
    if cram:
        ov["cram"] = str(_resolve(cram))
    if reference:
        ov["reference"] = str(_resolve(reference))
    if index:
        ov["index"] = str(_resolve(index))
    if min_mapq is not None:
        ov["min_mapping_quality"] = min_mapq
    if chroms:
        ov["chroms"] = tuple(chroms.split(","))
    if per_base:
        ov["per_base"] = True
    return ov


def _analysis_overrides(*, window=None, strata=None, per_base=False, outdir=None) -> dict:
    ov: dict = {}
    if window:
        ov["window"] = window
    if strata:
        ov["strata"] = dict(kv.split("=", 1) for kv in strata.split(",") if kv)
    if per_base:
        ov["per_base"] = True
    if outdir:
        ov["outdir"] = str(Path(outdir).expanduser().resolve())
    return ov


def _load_run(config, *, coverage=None, analysis=None) -> RunConfig:
    try:
        return RunConfig.load(config, coverage=coverage, analysis=analysis)
    except ValidationError as e:
        raise click.UsageError(
            "invalid or incomplete configuration -- provide --cram and --reference "
            f"(or a --config that sets them).\n{e}"
        )


def _echo_table(rows) -> None:
    click.echo("\t".join(TSV_COLUMNS))
    for r in rows:
        cells = r.as_row()
        click.echo("\t".join(str(cells[c]) for c in TSV_COLUMNS))


# Shared input options (composed onto coverage + analyze).
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


@click.group()
@click.version_option(package_name="parklab-chromcov", prog_name="chromcov")
def main() -> None:
    """Per-chromosome average coverage from a CRAM, with QC extensions."""


@main.command()
@input_options
@click.option("--chroms", default=None, help="comma-separated contig subset (default: all)")
@click.option("--output", type=click.Path(dir_okay=False), default=None,
              help="write the table to FILE (default: stdout)")
@click.option("--write", is_flag=True, help="also archive under runs/<name>/")
@click.option("--runs-dir", "runs_dir", type=click.Path(file_okay=False), default=None,
              help="archive directory (default: ./runs)")
@click.option("--run-name", "run_name", type=click.Choice(["slug", "hash"]), default="slug")
def coverage(cram, reference, index, min_mapq, config, chroms,
             output, write, runs_dir, run_name) -> None:
    """Per-chromosome mean coverage table (the deliverable)."""
    run = _load_run(config, coverage=_coverage_overrides(
        cram=cram, reference=reference, index=index, min_mapq=min_mapq, chroms=chroms))
    cfg = run.coverage
    rows = dispatch.run_coverage(cfg)

    if output:
        out = Path(output).expanduser().resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(rows, out)
        click.echo(f"wrote {out}", err=True)
    else:
        _echo_table(rows)

    if write:
        store = RunStore(_resolve(runs_dir) or Path("runs"))
        run_dir = store.write_run(rows, cfg, name_style=run_name)
        click.echo(f"wrote {run_dir}/coverage.tsv (+ .provenance.json)", err=True)


def _write_run_sidecar(path, run, config_file, store, analysis, akey, chrom_list) -> None:
    """run.json for a Level-2 analysis run. Self-describing + re-runnable: it
    embeds the *resolved* RunConfig (so the run reproduces from the sidecar alone)
    and points back at the coverage-key (Level-1 tracks) and the source config
    file, if any."""
    baseline, source = analysis.baseline()
    record = {
        "schema": "chromcov.analysis-run/2",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "tool": provenance.tool_version(),
        "code": provenance.git_provenance(),
        "analysis_key": akey,
        "coverage_key": store.key,
        "perbase_dir": str(store.dir),
        "config_file": str(_resolve(config_file)) if config_file else None,
        "config": run.model_dump(mode="json"),         # the whole resolved run config
        "chroms_run": chrom_list or "all",
        "chromosomes": analysis.chroms,
        "baseline": {"value": round(baseline, 4), "source": source},
        "flagged": {c: fl for c, fl in analysis.flagged},
    }
    Path(path).write_text(json.dumps(record, indent=2, sort_keys=True))


@main.command()
@input_options
@click.option("--chroms", default=None, help="comma-separated contig subset (default: whole genome)")
@click.option("--window", type=int, default=None, help="windowed-track bin size (bp)")
@click.option("--strata", default="",
              help="callability strata as label=bed[,label=bed]; e.g. "
                   "easy=SMaHT_easy_hg38.bed.gz,difficult=...,extreme=...")
@click.option("--per-base", "per_base", is_flag=True,
              help="write per-base depth tracks (Level 1) so later runs reuse them")
@click.option("--outdir", type=click.Path(file_okay=False), default=None,
              help="output root (default: ./out). Tracks -> <outdir>/<coverage-key>/, "
                   "each analysis -> <outdir>/<coverage-key>/<analysis-key>/")
def analyze(cram, reference, index, min_mapq, config, chroms, window,
            strata, per_base, outdir) -> None:
    """Full QC suite (Level 2): stats, windows, strata, plots (+ the coverage
    table), reusing per-base tracks when present. Each run nests under the
    coverage dataset it derives from, so runs (e.g. stratified vs not) accumulate
    beside the tracks for comparison."""
    run = _load_run(
        config,
        coverage=_coverage_overrides(cram=cram, reference=reference, index=index, min_mapq=min_mapq),
        analysis=_analysis_overrides(window=window, strata=strata, per_base=per_base, outdir=outdir),
    )
    cfg, acfg = run.coverage, run.analysis

    store = PerBaseStore(acfg.outdir, cfg)
    strata_obj = Strata.from_arg(acfg.strata)
    analysis = CoverageAnalysis(cfg, acfg, strata_obj)

    chrom_list = chroms.split(",") if chroms else None
    report = analysis.run(chroms=chrom_list, store=store, write_tracks=acfg.per_base)
    click.echo(f"[preflight] {report['reference_check']['status']}", err=True)

    labels = strata_obj.labels()
    akey = analysis_key(store.key, acfg, labels)
    run_dir = store.dir / analysis_slug(acfg, labels, akey)   # nest under the coverage-key
    run_dir.mkdir(parents=True, exist_ok=True)
    written = analysis.write_outputs(run_dir)
    _write_run_sidecar(run_dir / "run.json", run, config, store, analysis, akey, chrom_list)

    for line in analysis.summary_lines():
        click.echo(line, err=True)
    if acfg.per_base:
        click.echo(f"per-base tracks (Level 1): {store.dir}", err=True)
    click.echo(f"analysis run (Level 2): {run_dir}", err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)


@main.command()
@input_options
@click.option("--chrom", required=True, help="the single chromosome to compute a track for")
@click.option("--outdir", type=click.Path(file_okay=False), default="out",
              help="output root (default: out). Track -> <outdir>/<coverage-key>/")
def perbase(cram, reference, index, min_mapq, config, chrom, outdir) -> None:
    """Compute + store ONE chromosome's per-base depth track (Level 1).

    The Snakemake scatter unit; a no-op if the track already exists. Finalize the
    coverage.json sidecar by running `analyze --per-base` (the gather step)."""
    run = _load_run(config, coverage=_coverage_overrides(
        cram=cram, reference=reference, index=index, min_mapq=min_mapq, per_base=True))
    cfg = run.coverage
    store = PerBaseStore(Path(outdir).expanduser().resolve(), cfg)
    if store.has(chrom):
        click.echo(f"track exists (skip): {store.track_path(chrom)}", err=True)
        return
    summary = build_track(cfg, chrom, store)
    click.echo(f"wrote {store.track_path(chrom)}  "
               f"(mean {summary['mean']}x over {summary['length']:,} bp)", err=True)


@main.command()
@click.option("--runs-dir", "runs_dir", type=click.Path(file_okay=False), default=None,
              help="archive directory (default: ./runs)")
def collate(runs_dir) -> None:
    """Compare archived runs (wide chrom x run table)."""
    store = RunStore(_resolve(runs_dir) or Path("runs"))
    long_rows = store.collate()
    if not long_rows:
        click.echo(f"no runs under {store.runs_dir} yet (use `coverage --write`).")
        return
    run_ids, table = store.pivot_mean(long_rows)
    click.echo("chrom\t" + "\t".join(run_ids))
    for chrom, per_run in table.items():
        click.echo(chrom + "\t" + "\t".join(str(per_run.get(rid, "")) for rid in run_ids))
    click.echo("\n# run_id -> params")
    seen: dict = {}
    for r in long_rows:
        seen.setdefault(r["run_id"], {k: r[k] for k in ("min_mapping_quality",
                                                        "exclude_flags") if k in r})
    for rid, params in seen.items():
        click.echo(f"# {rid}: {params}")


@main.group(name="fetch")
def fetch_group() -> None:
    """Download inputs a clean clone needs (callability strata, ...)."""


@fetch_group.command(name="strata")
@click.option("--dest", type=click.Path(file_okay=False), default="data",
              help="directory to download into (default: ./data)")
@click.option("--force", is_flag=True, help="re-download even if the files exist")
def fetch_strata_cmd(dest, force) -> None:
    """Download the Park Lab SMaHT easy/difficult/extreme hg38 BEDs."""
    paths = fetch.fetch_strata(dest, force=force)
    click.echo("\nuse them with:\n  chromcov analyze ... --strata " + fetch.strata_arg(paths))


if __name__ == "__main__":
    main()
