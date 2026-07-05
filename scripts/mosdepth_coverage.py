#!/usr/bin/env python
"""
mosdepth_coverage.py — optional add-on: run mosdepth and emit the SAME outputs
chromcov does, so it's a drop-in cross-check rather than a second backend baked
into the tool.

    uv run python scripts/mosdepth_coverage.py --cram data/COLO829T_TEST.cram \
        --reference data/…_GRCh38_….fa --min-mapq 0 --outdir out/mosdepth

It shells out to the `mosdepth` binary on PATH (install: `mamba install -c
bioconda mosdepth`, or use the biocontainer), converts its
`.mosdepth.summary.txt` into a chromcov-format `coverage.tsv`
(chrom/length/bases/mean_coverage), and — with --per-base — splits mosdepth's
per-base BED into one `<chrom>.per-base.bedgraph.gz` per contig, matching
chromcov's track naming. Pass --compare <chromcov coverage.tsv> to print a
per-chromosome delta.

Why the numbers should agree: both count aligned bases over M/=/X CIGAR ops with
no overlap correction, and `mean = bases / length`. The default --exclude-flags
matches chromcov's pinned union mask (chromcov.config.DEFAULT_EXCLUDE), so the
two agree out of the box; mosdepth has no -G / exclude-all-flags equivalent, so
that knob is intentionally unsupported here.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

# The add-on borrows the package's output format so the results are byte-comparable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT))
from chromcov.config import DEFAULT_EXCLUDE, to_mask   # noqa: E402
from chromcov.result import ChromCoverage, write_tsv    # noqa: E402


def build_argv(cram, reference, prefix, *, min_mapq, exclude_flags, include_flags,
               threads, per_base, chrom) -> list[str]:
    argv = ["mosdepth"]
    if not per_base:
        argv.append("--no-per-base")
    argv += ["--threads", str(threads), "--fasta", str(reference),
             "--mapq", str(min_mapq), "--flag", str(exclude_flags)]
    if chrom:
        argv += ["--chrom", chrom]
    if include_flags:
        argv += ["--include-flag", str(include_flags)]
    argv += [str(prefix), str(cram)]
    return argv


def parse_summary(summary_path: Path) -> list[ChromCoverage]:
    """<prefix>.mosdepth.summary.txt -> chromcov rows. Skips the `total` row and
    any `_region` rows; mean is recomputed as bases/length to match chromcov."""
    rows: list[ChromCoverage] = []
    with open(summary_path) as fh:
        for rec in csv.DictReader(fh, delimiter="\t"):
            chrom = rec["chrom"]
            if chrom == "total" or chrom.endswith("_region"):
                continue
            rows.append(ChromCoverage(chrom=chrom, length=int(rec["length"]), bases=int(rec["bases"])))
    return rows


def split_per_base(perbase_bed: Path, outdir: Path) -> list[Path]:
    """Split mosdepth's combined PREFIX.per-base.bed.gz into one
    <chrom>.per-base.bedgraph.gz per contig (chromcov's per-chrom track naming)."""
    written: list[Path] = []
    current = None
    fh_out = None
    try:
        with gzip.open(perbase_bed, "rt") as fh:
            for line in fh:
                chrom = line.split("\t", 1)[0]
                if chrom != current:
                    if fh_out:
                        fh_out.close()
                    current = chrom
                    path = outdir / f"{chrom}.per-base.bedgraph.gz"
                    fh_out = gzip.open(path, "wt")
                    written.append(path)
                fh_out.write(line)
    finally:
        if fh_out:
            fh_out.close()
    return written


def compare(rows: list[ChromCoverage], chromcov_tsv: Path) -> None:
    ours = {}
    with open(chromcov_tsv) as fh:
        for rec in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            ours[rec["chrom"]] = float(rec["mean_coverage"])
    print(f"\n# comparison vs {chromcov_tsv}")
    print("chrom\tchromcov\tmosdepth\tdelta")
    for r in rows:
        c = ours.get(r.chrom)
        if c is None:
            continue
        print(f"{r.chrom}\t{c:.2f}\t{r.mean:.2f}\t{r.mean - c:+.2f}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Run mosdepth and emit chromcov-format outputs.")
    ap.add_argument("--cram", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--min-mapq", type=int, default=0)
    ap.add_argument("--exclude-flags", default=str(DEFAULT_EXCLUDE),
                    help="int mask or comma-separated flag names (default matches chromcov)")
    ap.add_argument("--include-flags", default="0", help="int mask or comma-separated flag names")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chrom", default=None, help="restrict to a single contig (mosdepth -c)")
    ap.add_argument("--per-base", action="store_true", help="also split per-base tracks")
    ap.add_argument("--outdir", default="out/mosdepth")
    ap.add_argument("--compare", default=None, help="a chromcov coverage.tsv to diff against")
    args = ap.parse_args()

    if shutil.which("mosdepth") is None:
        sys.exit("mosdepth not found on PATH (mamba install -c bioconda mosdepth, or the biocontainer).")

    exclude = to_mask([f.strip() for f in args.exclude_flags.split(",")]
                      if not args.exclude_flags.lstrip("-").isdigit() else int(args.exclude_flags))
    include = to_mask([f.strip() for f in args.include_flags.split(",")]
                      if not args.include_flags.lstrip("-").isdigit() else int(args.include_flags))

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / "cov"
        argv = build_argv(args.cram, args.reference, prefix, min_mapq=args.min_mapq,
                          exclude_flags=exclude, include_flags=include, threads=args.threads,
                          per_base=args.per_base, chrom=args.chrom)
        print("[mosdepth] " + " ".join(argv), file=sys.stderr)
        subprocess.run(argv, check=True)

        rows = parse_summary(prefix.with_suffix(".mosdepth.summary.txt"))
        coverage_tsv = outdir / "coverage.tsv"
        write_tsv(rows, coverage_tsv)
        print(f"wrote {coverage_tsv} ({len(rows)} chromosomes)", file=sys.stderr)

        if args.per_base:
            tracks = split_per_base(prefix.with_suffix(".per-base.bed.gz"), outdir)
            print(f"wrote {len(tracks)} per-base tracks to {outdir}/", file=sys.stderr)

    if args.compare:
        compare(rows, Path(args.compare))


if __name__ == "__main__":
    main()
