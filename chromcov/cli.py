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

from pathlib import Path

import click

from . import dispatch
from .config import AnalysisConfig, CoverageConfig
from .output import RunStore
from .pipeline import CoverageAnalysis
from .result import TSV_COLUMNS, write_tsv
from .strata import Strata


def _resolve(p: str | None) -> Path | None:
    return Path(p).expanduser().resolve() if p else None


def _build_coverage_config(cram, reference, index, min_mapq, config,
                           backend=None, chroms=None, per_base=False) -> CoverageConfig:
    """Build a validated CoverageConfig from CLI options, optionally layered over
    a --config YAML (explicit options win)."""
    overrides: dict = {}
    if cram:
        overrides["cram"] = str(_resolve(cram))
    if reference:
        overrides["reference"] = str(_resolve(reference))
    if index:
        overrides["index"] = str(_resolve(index))
    if backend:
        overrides["backend"] = backend
    if min_mapq is not None:
        overrides["min_mapping_quality"] = min_mapq
    if chroms:
        overrides["chroms"] = tuple(chroms.split(","))
    if per_base:
        overrides["per_base"] = True

    if config:
        base = CoverageConfig.from_yaml(config).model_dump()
        base.update(overrides)
        return CoverageConfig.model_validate(base)
    if "cram" not in overrides or "reference" not in overrides:
        raise click.UsageError("--cram and --reference are required (or pass --config).")
    return CoverageConfig.model_validate(overrides)


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
@click.option("--backend", type=click.Choice(["native", "mosdepth"]), default=None)
@click.option("--chroms", default=None, help="comma-separated contig subset (default: all)")
@click.option("--output", type=click.Path(dir_okay=False), default=None,
              help="write the table to FILE (default: stdout)")
@click.option("--write", is_flag=True, help="also archive under runs/<name>/")
@click.option("--runs-dir", "runs_dir", type=click.Path(file_okay=False), default=None,
              help="archive directory (default: ./runs)")
@click.option("--run-name", "run_name", type=click.Choice(["slug", "hash"]), default="slug")
def coverage(cram, reference, index, min_mapq, config, backend, chroms,
             output, write, runs_dir, run_name) -> None:
    """Per-chromosome mean coverage table (the deliverable)."""
    cfg = _build_coverage_config(cram, reference, index, min_mapq, config,
                                 backend=backend, chroms=chroms)
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


@main.command()
@input_options
@click.option("--chroms", default=None, help="comma-separated contig subset (default: whole genome)")
@click.option("--window", type=int, default=None, help="windowed-track bin size (bp)")
@click.option("--strata", default="",
              help="callability strata as label=bed[,label=bed]; e.g. "
                   "easy=SMaHT_easy_hg38.bed.gz,difficult=...,extreme=...")
@click.option("--per-base", "per_base", is_flag=True,
              help="also write the RLE per-base bedgraph.gz (bulky on WGS)")
@click.option("--outdir", type=click.Path(file_okay=False), default=None,
              help="output directory (default: ./out)")
def analyze(cram, reference, index, min_mapq, config, chroms, window,
            strata, per_base, outdir) -> None:
    """Full QC suite: stats, windows, strata, plots (+ the coverage table)."""
    cfg = _build_coverage_config(cram, reference, index, min_mapq, config)

    acfg = AnalysisConfig.from_yaml(config) if config else AnalysisConfig()
    if window:
        acfg.window = window
    if outdir:
        acfg.outdir = Path(outdir).expanduser().resolve()
    if per_base:
        acfg.per_base = True

    strata_obj = Strata.from_arg(strata) if strata else Strata.from_arg(acfg.strata)

    analysis = CoverageAnalysis(cfg, acfg, strata_obj)
    out = acfg.outdir
    out.mkdir(parents=True, exist_ok=True)

    chrom_list = chroms.split(",") if chroms else None
    per_base_path = out / "coverage.perbase.bedgraph.gz" if acfg.per_base else None

    report = analysis.run(chroms=chrom_list, per_base_path=per_base_path)
    click.echo(f"[preflight] ok: sorted, indexed, reference {report['reference_check']['status']}",
               err=True)

    written = analysis.write_outputs(out, per_base=acfg.per_base)
    for line in analysis.summary_lines():
        click.echo(line, err=True)
    click.echo(f"wrote {len(written)} outputs to {out}/", err=True)
    for name, path in written.items():
        click.echo(f"  {name}: {path.name}", err=True)


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
        seen.setdefault(r["run_id"], {k: r[k] for k in ("backend", "min_mapping_quality",
                                                         "exclude_flags") if k in r})
    for rid, params in seen.items():
        click.echo(f"# {rid}: {params}")


if __name__ == "__main__":
    main()
