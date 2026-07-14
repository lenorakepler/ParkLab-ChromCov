#!/usr/bin/env python
"""
mosdepth_coverage.py — optional add-on: run mosdepth and emit the SAME outputs
chromcov does, so it's a drop-in cross-check rather than a second backend baked
into the tool.

    # cross-check an existing chromcov --full run (matches its settings + baseline):
    uv run python scripts/mosdepth_coverage.py --run out

    # or drive it directly:
    uv run python scripts/mosdepth_coverage.py --cram data/COLO829T_TEST.cram \
        --reference data/…_GRCh38_….fa --min-mapq 0 --outdir out/mosdepth

It shells out to the `mosdepth` binary on PATH (install: `mamba install -c
bioconda mosdepth`, or use the biocontainer), converts its
`.mosdepth.summary.txt` into a chromcov-format `coverage.tsv`
(chrom/length/bases/mean_coverage), and — with --per-base — splits mosdepth's
per-base BED into one `<chrom>.per-base.bedgraph.gz` per contig, matching
chromcov's track naming.

--run <chromcov --full output dir or its run.json> makes the comparison
genuinely apples-to-apples: it reads the run's embedded Config and pulls the
parameters that determine the numbers — reference/CRAM, min mapping quality,
include/exclude flag masks, window size, ploidy, scatter CN cap — so mosdepth is
invoked with the same filtering, and diffs against that run's coverage.tsv by
default. Any explicit CLI flag still overrides the value from the run.

Plots are emitted by default (--no-plots to skip), reusing chromcov's own
plotting code (chromcov.present.plots) so the figures line up with a chromcov
run: the per-chromosome bar chart always, and the windowed copy-number scatter +
ideogram whenever a window size is known (from the run's `window`, or --by),
which makes mosdepth also emit a `.regions.bed.gz` of per-window means. With
--run the copy-number baseline is chromcov's *own* computed baseline value (read
from run.json), so the CN axis is identical; without a run it falls back to the
median of autosomal means (of windows if available, else whole chromosomes) — a
coarse stand-in for chromcov's easy-autosomal-median, which needs a per-base
histogram mosdepth's summary doesn't carry. mosdepth has no callability strata,
so scatter points render in one color rather than easy/difficult/extreme tiers.

How the numbers relate: both tally aligned bases over each read's gapless M/=/X
blocks and take `mean = bases / length`, under the same flag mask — but they are
NOT identical, and that is the point of putting the plots side by side. chromcov
reads blocks from pysam `AlignedSegment.get_blocks()`, which is per-read and does
not deduplicate overlapping mate pairs, whereas mosdepth (its default mode, used
here) decrements the overlap of a fragment's two mates so it counts once. So
chromcov reads slightly HIGH wherever mates overlap; the bar/scatter comparison
makes that bias visible (chromcov ≥ mosdepth) until chromcov grows its own
mate-overlap correction. mosdepth has no -G / exclude-all-flags equivalent, so
that knob is unsupported here (a --run that set it is warned about).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

# The add-on borrows the package's output format so the results are byte-comparable.
_REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_REPO_ROOT))
from chromcov.filtering import DEFAULT_EXCLUDE, to_mask   # noqa: E402
from chromcov.policy import is_autosome   # noqa: E402
from chromcov.present.frames import ChromCoverage, coverage_frame, write_table    # noqa: E402


def resolve_mask(value) -> int:
    """Accept a flag mask as an int (as run.json stores it) or a CLI string that is
    either an int literal or comma-separated flag names; return the int mask."""
    if isinstance(value, int):
        return to_mask(value)
    s = str(value)
    return to_mask(int(s) if s.lstrip("-").isdigit()
                   else [f.strip() for f in s.split(",")])


def load_run(path: str | Path) -> dict:
    """Read a chromcov --full run (a run.json, or the dir holding it) and pull the
    settings that determine the coverage numbers, so this cross-check runs mosdepth
    the same way and diffs against the same baseline. Returns a dict of seeds the
    CLI layers explicit overrides on top of."""
    p = Path(path).expanduser().resolve()
    sidecar = p / "run.json" if p.is_dir() else p
    rec = json.loads(sidecar.read_text())
    cfg = rec.get("config", {})
    chroms = cfg.get("chroms")
    return {
        "cram": cfg.get("cram"),
        "reference": cfg.get("reference"),
        "min_mapq": cfg.get("min_mapping_quality", 0),
        "exclude_flags": cfg.get("exclude_flags", DEFAULT_EXCLUDE),
        "include_flags": cfg.get("include_flags", 0),
        "exclude_all_flags": cfg.get("exclude_all_flags", 0),
        "window": cfg.get("window"),
        "ploidy": cfg.get("ploidy", 2),
        "cap_cn": cfg.get("scatter_cap_cn", 6.0),
        # chromcov's own computed diploid reference -> identical CN axis on the plots.
        "baseline": (rec.get("baseline") or {}).get("value"),
        # mosdepth --chrom restricts to a single contig only; a multi-contig subset
        # is left to run over all contigs (compare() only shows the overlap).
        "chrom": chroms[0] if chroms and len(chroms) == 1 else None,
        "coverage_tsv": str(sidecar.parent / "coverage.tsv"),
    }


def build_argv(cram, reference, prefix, *, min_mapq, exclude_flags, include_flags,
               threads, per_base, chrom, by=None) -> list[str]:
    argv = ["mosdepth"]
    if not per_base:
        argv.append("--no-per-base")
    if by:
        argv += ["--by", str(by)]   # int window size or a BED path -> PREFIX.regions.bed.gz
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


def parse_regions(regions_bed: Path) -> list[dict]:
    """mosdepth's PREFIX.regions.bed.gz (from --by) -> window rows for the scatter:
    {chrom, start, end, mean}. No strata/easy_frac (mosdepth has none), so the
    scatter draws every window in a single color."""
    windows: list[dict] = []
    with gzip.open(regions_bed, "rt") as fh:
        for line in fh:
            chrom, start, end, mean = line.rstrip("\n").split("\t")[:4]
            windows.append({"chrom": chrom, "start": int(start), "end": int(end),
                            "mean": float(mean)})
    return windows


def autosomal_baseline(pairs) -> float:
    """Median over autosomes of an iterable of (chrom, value) pairs -- a coarse
    stand-in for chromcov's easy-autosomal-median baseline (mosdepth's summary
    carries no per-base histogram to take an easy-masked median from). Fed
    per-chromosome means for the bar, or per-window means for the scatter."""
    auto = [v for c, v in pairs if is_autosome(c)]
    return statistics.median(auto) if auto else 0.0


def make_plots(rows: list[ChromCoverage], windows: list[dict] | None, outdir: Path,
               *, baseline: float | None = None, ploidy: int = 2,
               cap_cn: float = 6.0) -> list[Path]:
    """Reuse chromcov's plotting on the mosdepth numbers so the figures line up
    with a chromcov run. Bar always; scatter only when a window size produced
    `windows`. `baseline` is chromcov's own computed value (from --run) so the CN
    axis is identical; when None, fall back to the autosomal-median stand-in at the
    granularity each plot needs."""
    from chromcov.present import plots   # lazy: pulls in matplotlib/plotly/kaleido

    written: list[Path] = []
    chrom_means = {r.chrom: r.mean for r in rows}
    bar_baseline = baseline or autosomal_baseline((r.chrom, r.mean) for r in rows)

    bar = outdir / "coverage.bar.png"
    plots.bar_by_chromosome(chrom_means, bar, baseline=bar_baseline or None)
    written.append(bar)

    if windows:
        scatter_baseline = baseline or autosomal_baseline((w["chrom"], w["mean"]) for w in windows)
        scatter = outdir / "coverage.scatter.png"
        out = plots.scatter_windows(windows, scatter, baseline=scatter_baseline,
                                    ploidy=ploidy, cap_cn=cap_cn)
        written += [out["png"], out["html"]]
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
    ap.add_argument("--run", default=None,
                    help="a chromcov --full output dir (or its run.json): pull cram/reference/"
                         "filters/window/ploidy/baseline so the comparison matches, and diff "
                         "against its coverage.tsv by default")
    ap.add_argument("--cram", default=None, help="the CRAM to measure (overrides the run's)")
    ap.add_argument("--reference", default=None, help="reference FASTA (overrides the run's)")
    ap.add_argument("--min-mapq", type=int, default=None, help="-Q override (default: the run's, else 0)")
    ap.add_argument("--exclude-flags", default=None,
                    help="int mask or comma-separated flag names (default: the run's, else chromcov's)")
    ap.add_argument("--include-flags", default=None, help="int mask or comma-separated flag names")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--chrom", default=None, help="restrict to a single contig (mosdepth -c)")
    ap.add_argument("--by", default=None,
                    help="windowed-mean bin size for the scatter (mosdepth --by): int bp (e.g. "
                         "10000) or a BED path; default: the run's `window`")
    ap.add_argument("--per-base", action="store_true", help="also split per-base tracks")
    ap.add_argument("--no-plots", action="store_true", help="skip the bar/scatter plots")
    ap.add_argument("--outdir", default="out/mosdepth")
    ap.add_argument("--compare", default=None,
                    help="a chromcov coverage.tsv to diff against (default: the run's)")
    args = ap.parse_args()

    if shutil.which("mosdepth") is None:
        sys.exit("mosdepth not found on PATH (mamba install -c bioconda mosdepth, or the biocontainer).")

    # Config precedence: an explicit CLI flag wins, else the --run value, else the default.
    run = load_run(args.run) if args.run else {}
    cram = args.cram or run.get("cram")
    reference = args.reference or run.get("reference")
    if not (cram and reference):
        sys.exit("need --cram and --reference (or a --run that provides them)")

    min_mapq = args.min_mapq if args.min_mapq is not None else run.get("min_mapq", 0)
    exclude = resolve_mask(args.exclude_flags if args.exclude_flags is not None
                           else run.get("exclude_flags", DEFAULT_EXCLUDE))
    include = resolve_mask(args.include_flags if args.include_flags is not None
                           else run.get("include_flags", 0))
    by = args.by if args.by is not None else run.get("window")
    chrom = args.chrom or run.get("chrom")
    compare_tsv = args.compare or run.get("coverage_tsv")

    if run.get("exclude_all_flags"):
        print(f"[warn] the run set exclude_all_flags (-G {run['exclude_all_flags']}); mosdepth "
              "has no -G equivalent, so that filter is NOT reproduced here.", file=sys.stderr)

    want_plots = not args.no_plots
    if want_plots and not by:
        print("[note] no window size (--by, or a run's `window`): bar chart only, no scatter.",
              file=sys.stderr)

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as td:
        prefix = Path(td) / "cov"
        argv = build_argv(cram, reference, prefix, min_mapq=min_mapq,
                          exclude_flags=exclude, include_flags=include, threads=args.threads,
                          per_base=args.per_base, chrom=chrom, by=by)
        print("[mosdepth] " + " ".join(argv), file=sys.stderr)
        subprocess.run(argv, check=True)

        rows = parse_summary(prefix.with_suffix(".mosdepth.summary.txt"))
        coverage_tsv = outdir / "coverage.tsv"
        write_table(coverage_frame(rows), coverage_tsv)
        print(f"wrote {coverage_tsv} ({len(rows)} chromosomes)", file=sys.stderr)

        if args.per_base:
            tracks = split_per_base(prefix.with_suffix(".per-base.bed.gz"), outdir)
            print(f"wrote {len(tracks)} per-base tracks to {outdir}/", file=sys.stderr)

        if want_plots:
            windows = parse_regions(prefix.with_suffix(".regions.bed.gz")) if by else None
            figs = make_plots(rows, windows, outdir, baseline=run.get("baseline"),
                              ploidy=run.get("ploidy", 2), cap_cn=run.get("cap_cn", 6.0))
            print(f"wrote {len(figs)} plot file(s): {', '.join(p.name for p in figs)}",
                  file=sys.stderr)

    if compare_tsv and Path(compare_tsv).exists():
        compare(rows, Path(compare_tsv))
    elif args.compare:
        sys.exit(f"--compare target not found: {args.compare}")


if __name__ == "__main__":
    main()
