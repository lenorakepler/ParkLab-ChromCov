"""
End-to-end smoke run of the dispatcher against data/, with run archival + collate.

    # one run, archived under runs/<hash>/
    uv run python dev/dispatch-sketch/run_example.py --backend native --chroms chr21 --write

    # a second run with a different filter -> a different runs/<hash>/
    uv run python dev/dispatch-sketch/run_example.py --backend native --chroms chr21 --min-mapq 20 --write

    # compare everything under runs/ (wide chrom x run table)
    uv run python dev/dispatch-sketch/run_example.py --collate

`--chroms` restricts the native run to a few contigs so this stays quick on the
17 GB CRAM (the full genome loop is the same code path, just slower). The
mosdepth backend needs the binary on PATH.
"""
from __future__ import annotations
import argparse
from pathlib import Path

from config import CoverageConfig
from result import TSV_COLUMNS
import dispatch
import output

_REPO_ROOT = Path(__file__).resolve().parents[2]
_DATA = _REPO_ROOT / "data"
_RUNS = _REPO_ROOT / "runs"


def _print_collation() -> None:
    long_rows = output.collate(_RUNS)
    if not long_rows:
        print(f"no runs under {_RUNS} yet (use --write).")
        return
    run_ids, table = output.pivot_mean(long_rows)
    print("chrom\t" + "\t".join(run_ids))
    for chrom, per_run in table.items():
        print(chrom + "\t" + "\t".join(str(per_run.get(rid, "")) for rid in run_ids))
    print("\n# run_id -> params")
    seen = {}
    for r in long_rows:
        seen.setdefault(r["run_id"], {k: r[k] for k in output.RUN_PARAM_FIELDS if k in r})
    for rid, params in seen.items():
        print(f"# {rid}: {params}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", default="native", choices=["native", "mosdepth"])
    ap.add_argument("--chroms", default="", help="comma-separated subset (native only)")
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--write", action="store_true", help="archive under runs/<name>/")
    ap.add_argument("--run-name", default="slug", choices=["slug", "hash"],
                    help="run dir naming: slug (readable + hash suffix) or bare hash")
    ap.add_argument("--collate", action="store_true", help="compare all runs under runs/ and exit")
    args = ap.parse_args()

    if args.collate:
        _print_collation()
        return

    config = CoverageConfig(
        cram=_DATA / "COLO829T_TEST.cram",
        reference=_DATA / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa",
        backend=args.backend,
        min_mapping_quality=args.min_mapq,
        chroms=tuple(args.chroms.split(",")) if args.chroms else None,
    )

    rows = dispatch.run_coverage(config)

    print("\t".join(TSV_COLUMNS))
    for r in rows:
        cells = r.as_row()
        print("\t".join(str(cells[c]) for c in TSV_COLUMNS))

    if args.write:
        run_dir = output.write_run(rows, config, runs_dir=_RUNS, name_style=args.run_name)
        print(f"\nwrote {run_dir}/coverage.tsv (+ .provenance.json)")


if __name__ == "__main__":
    main()
