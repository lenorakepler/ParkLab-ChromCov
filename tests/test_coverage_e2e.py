"""
End-to-end coverage correctness against a synthetic CRAM with hand-computed depth.

This is the behavioral safety net for the Plan A refactor (dev/PLAN-A-DESIGN.md
step 1): it pins the *coverage numbers* on known input, so the mechanical moves
that follow (perbase demotion, map/reduce split, parallelism, resume, the polars
swap) can be verified to leave the reported values unchanged.

The fixture builds a tiny 2-contig reference and a CRAM whose reads deliberately
exercise the cases that are easy to get wrong:

  * a deletion (CIGAR D) and a ref-skip (CIGAR N) -- the skipped reference
    positions must contribute ZERO depth (get_blocks excludes D/N),
  * a soft-clip (CIGAR S) -- clipped bases must not count and must not shift the
    aligned block,
  * reads flagged secondary / supplementary -- dropped by the default filter; and a
    duplicate read -- KEPT by the default (config.DEFAULT_EXCLUDE), dropped only by
    an explicit drop-all filter.

chr1 (length 30) kept reads and their aligned blocks (0-based, half-open):
    A  pos 0  10M           -> [0,10)
    D  pos 2  3S6M          -> [2,8)            (soft-clip dropped)
    B  pos 5  5M3D5M        -> [5,10) [13,18)   (del 10,11,12 = 0 depth)
    C  pos 20 4M2N4M        -> [20,24) [26,30)  (skip 24,25 = 0 depth)
  dropped by the default filter: E(dup), F(secondary), G(supplementary), all 10M @0

Hand-computed per-base depth over chr1 positions 0..29 -> total 34 aligned bases.
chr2 (length 10): one 5M read @0 -> depth 1 over [0,5), total 5 aligned bases.
"""
from __future__ import annotations

import numpy as np
import pysam
import pytest

from click.testing import CliRunner

from chromcov import perbase
from chromcov.calc_cov import calc_cov
from chromcov.cli import main
from chromcov.config import Config, ReadFilter
from chromcov.coverage import run_coverage
from chromcov.qc_report import QCReport, compute_chrom

CHR1_LEN = 30
CHR2_LEN = 10

# Explicitly drop all three artifact classes, so the depth-math assertions are
# independent of the default policy (which keeps duplicates -- see the default
# test below). run_coverage/QCReport still use config.DEFAULT_EXCLUDE.
DROP_ALL = ["unmapped", "secondary", "duplicate", "supplementary"]

# The single source of truth the whole test asserts against.
EXPECTED_CHR1_DEPTH = np.array(
    [1, 1, 2, 2, 2, 3, 3, 3, 2, 2,   # 0..9
     0, 0, 0, 1, 1, 1, 1, 1, 0, 0,   # 10..19  (10,11,12 = deletion; 18,19 uncovered)
     1, 1, 1, 1, 0, 0, 1, 1, 1, 1],  # 20..29  (24,25 = ref skip)
    dtype=np.int32,
)
EXPECTED_CHR1_BASES = int(EXPECTED_CHR1_DEPTH.sum())   # 34 (A+B+C+D, drop-all filter)
EXPECTED_CHR2_BASES = 5
DUP_BASES = 10                                              # E_dup (10M @0): a duplicate read
EXPECTED_CHR1_BASES_DEFAULT = EXPECTED_CHR1_BASES + DUP_BASES   # 44: the default keeps duplicates


def _ref_seq(length: int) -> str:
    return ("ACGT" * ((length // 4) + 1))[:length]


def _seg(header, name, flag, tid, pos, cigar, seq):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.flag = flag
    a.reference_id = tid
    a.reference_start = pos
    a.mapping_quality = 60
    a.cigartuples = cigar
    a.query_sequence = seq
    a.next_reference_id = -1
    a.next_reference_start = -1
    a.template_length = 0
    return a


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    """Build a sorted, indexed synthetic CRAM + reference; return a Config."""
    d = tmp_path_factory.mktemp("synth")
    ref_path = d / "ref.fa"
    chr1 = _ref_seq(CHR1_LEN)
    chr2 = _ref_seq(CHR2_LEN)
    ref_path.write_text(f">chr1\n{chr1}\n>chr2\n{chr2}\n")
    pysam.faidx(str(ref_path))

    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": CHR1_LEN}, {"SN": "chr2", "LN": CHR2_LEN}],
    })

    M, D, N, S = 0, 2, 3, 4  # CIGAR op codes (I=1 unused here)
    DUP, SEC, SUPP = 0x400, 0x100, 0x800
    # Written in coordinate-sorted order (chr1 by pos, then chr2).
    reads = [
        _seg(header, "A", 0, 0, 0, [(M, 10)], chr1[0:10]),
        _seg(header, "E_dup", DUP, 0, 0, [(M, 10)], chr1[0:10]),          # dropped
        _seg(header, "F_secondary", SEC, 0, 0, [(M, 10)], chr1[0:10]),    # dropped
        _seg(header, "G_supp", SUPP, 0, 0, [(M, 5)], chr1[0:5]),          # dropped
        _seg(header, "D", 0, 0, 2, [(S, 3), (M, 6)], "TTT" + chr1[2:8]),
        _seg(header, "B", 0, 0, 5, [(M, 5), (D, 3), (M, 5)], chr1[5:10] + chr1[13:18]),
        _seg(header, "C", 0, 0, 20, [(M, 4), (N, 2), (M, 4)], chr1[20:24] + chr1[26:30]),
        _seg(header, "H", 0, 1, 0, [(M, 5)], chr2[0:5]),
    ]

    cram_path = d / "synthetic.cram"
    with pysam.AlignmentFile(str(cram_path), "wc", header=header,
                             reference_filename=str(ref_path)) as out:
        for r in reads:
            out.write(r)
    pysam.index(str(cram_path))

    return Config(cram=cram_path, reference=ref_path)


def _open(cfg):
    return pysam.AlignmentFile(str(cfg.cram), "rc",
                               reference_filename=str(cfg.reference),
                               index_filename=str(cfg.index))


def test_per_base_depth_exact(synthetic):
    """The per-base vector matches the hand-computed depth exactly, including the
    zeros at the deletion and ref-skip -- the load-bearing correctness claim."""
    cram = _open(synthetic)
    try:
        base_depth, total, cov = calc_cov(cram, "chr1", ReadFilter(exclude_flags=DROP_ALL), per_base=True)
    finally:
        cram.close()
    assert np.array_equal(base_depth, EXPECTED_CHR1_DEPTH)
    assert int(total) == EXPECTED_CHR1_BASES
    assert cov == pytest.approx(EXPECTED_CHR1_BASES / CHR1_LEN)


def test_aggregate_matches_per_base(synthetic):
    """The mean-only path (per_base=False) must sum to the same total as the
    per-base path -- the two are the same number computed two ways."""
    cram = _open(synthetic)
    try:
        rf = ReadFilter(exclude_flags=DROP_ALL)
        _, agg_total, agg_cov = calc_cov(cram, "chr1", rf, per_base=False)
        _, pb_total, _ = calc_cov(cram, "chr1", rf, per_base=True)
    finally:
        cram.close()
    assert int(agg_total) == int(pb_total) == EXPECTED_CHR1_BASES
    assert agg_cov == pytest.approx(EXPECTED_CHR1_BASES / CHR1_LEN)


def test_clearing_excludes_counts_all_reads(synthetic):
    """With excludes cleared, every read is counted -- the filter, not fetch, is
    what drops secondary/supplementary by default."""
    cram = _open(synthetic)
    try:
        keep_none = ReadFilter(exclude_flags=0)
        _, total_all, _ = calc_cov(cram, "chr1", keep_none, per_base=False)
    finally:
        cram.close()
    # A/B/C/D (34) + E_dup(10) + F_secondary(10) + G_supplementary(5) = 59.
    assert int(total_all) == EXPECTED_CHR1_BASES + 25


def test_default_keeps_duplicates(synthetic):
    """The shipped default (config.DEFAULT_EXCLUDE) drops secondary/supplementary
    but KEEPS duplicates -- a deliberate choice for high-depth cancer data -- so it
    counts the E_dup read that an explicit drop-all filter removes."""
    cram = _open(synthetic)
    try:
        _, default_total, _ = calc_cov(cram, "chr1", ReadFilter(), per_base=False)
        _, dropall_total, _ = calc_cov(cram, "chr1", ReadFilter(exclude_flags=DROP_ALL), per_base=False)
    finally:
        cram.close()
    assert int(default_total) == EXPECTED_CHR1_BASES_DEFAULT   # 44: keeps the duplicate
    assert int(dropall_total) == EXPECTED_CHR1_BASES           # 34: duplicate removed
    assert int(default_total) - int(dropall_total) == DUP_BASES


def test_run_coverage_table(synthetic):
    """The whole-genome mean table (through preflight + contig selection) reports
    the hand-computed means for both contigs."""
    rows = run_coverage(synthetic)
    by_chrom = {r.chrom: r for r in rows}
    assert set(by_chrom) == {"chr1", "chr2"}

    c1 = by_chrom["chr1"]
    assert c1.length == CHR1_LEN
    assert c1.bases == EXPECTED_CHR1_BASES_DEFAULT   # default keeps the duplicate read
    assert c1.mean == pytest.approx(EXPECTED_CHR1_BASES_DEFAULT / CHR1_LEN)

    c2 = by_chrom["chr2"]
    assert c2.length == CHR2_LEN
    assert c2.bases == EXPECTED_CHR2_BASES
    assert c2.mean == pytest.approx(EXPECTED_CHR2_BASES / CHR2_LEN)


# --- --full path: compute (bedgraphs) -> gather (stats/plots) ----------------

def test_full_matches_mean_and_writes_outputs(synthetic, tmp_path):
    """The --full per-base pass reports the SAME mean as the mean-only path, plus
    stats, and writes the per-chrom bedgraphs + the analysis outputs."""
    bedgraph_dir = tmp_path / "perbase"
    rep = QCReport(synthetic)
    rep.run(bedgraph_dir=bedgraph_dir, jobs=1)

    rows = {r.chrom: r for r in rep.coverage_rows()}
    assert rows["chr1"].bases == EXPECTED_CHR1_BASES_DEFAULT
    assert rows["chr1"].mean == pytest.approx(EXPECTED_CHR1_BASES_DEFAULT / CHR1_LEN)
    assert rows["chr1"].stats is not None            # full carries per-base stats
    assert perbase.has_bedgraph(bedgraph_dir, "chr1")
    assert perbase.has_bedgraph(bedgraph_dir, "chr2")

    written = rep.write_outputs(tmp_path)
    assert (tmp_path / "coverage.tsv").exists()
    assert written["bar"].exists() and written["scatter"].exists()


def test_resume_reuses_existing_bedgraph(synthetic, tmp_path):
    """A contig already on disk is not recomputed (resume), and the run still
    produces the full genome."""
    bedgraph_dir = tmp_path / "perbase"
    compute_chrom(synthetic, "chr1", bedgraph_dir)          # pretend a prior partial run
    mtime = perbase.bedgraph_path(bedgraph_dir, "chr1").stat().st_mtime_ns

    QCReport(synthetic).run(bedgraph_dir=bedgraph_dir, jobs=1)
    assert perbase.bedgraph_path(bedgraph_dir, "chr1").stat().st_mtime_ns == mtime  # untouched
    assert perbase.has_bedgraph(bedgraph_dir, "chr2")       # the rest got computed


def test_gather_is_cram_free_and_reproduces(synthetic, tmp_path):
    """`plot`/gather reduces from the bedgraphs alone (lengths from the reference),
    reproducing the same per-chrom bases."""
    bedgraph_dir = tmp_path / "perbase"
    QCReport(synthetic).run(bedgraph_dir=bedgraph_dir, jobs=1)

    regathered = QCReport(synthetic)
    regathered.gather(bedgraph_dir)                          # no CRAM used
    rows = {r.chrom: r.bases for r in regathered.coverage_rows()}
    assert rows == {"chr1": EXPECTED_CHR1_BASES_DEFAULT, "chr2": EXPECTED_CHR2_BASES}


def test_jobs_parallel_matches_serial(synthetic, tmp_path):
    serial = QCReport(synthetic)
    serial.run(bedgraph_dir=tmp_path / "s", jobs=1)
    parallel = QCReport(synthetic)
    parallel.run(bedgraph_dir=tmp_path / "p", jobs=2)
    assert {r.chrom: r.bases for r in serial.coverage_rows()} == \
           {r.chrom: r.bases for r in parallel.coverage_rows()}


# --- CLI surface -------------------------------------------------------------

def test_cli_mean_default_to_stdout(synthetic):
    res = CliRunner().invoke(main, ["coverage", "--cram", str(synthetic.cram),
                                    "--reference", str(synthetic.reference)])
    assert res.exit_code == 0, res.output
    assert "mean_coverage" in res.output and "chr1" in res.output


def test_cli_full_then_plot(synthetic, tmp_path):
    runner = CliRunner()
    res = runner.invoke(main, ["coverage", "--cram", str(synthetic.cram),
                               "--reference", str(synthetic.reference),
                               "--full", "--outdir", str(tmp_path)])
    assert res.exit_code == 0, res.output
    assert (tmp_path / "coverage.tsv").exists()
    assert (tmp_path / "run.json").exists()
    assert perbase.has_bedgraph(tmp_path / "perbase", "chr1")

    res = runner.invoke(main, ["plot", "--outdir", str(tmp_path)])   # CRAM-free re-graph
    assert res.exit_code == 0, res.output
    assert (tmp_path / "coverage.bar.png").exists()


def test_windows_tagged_with_dominant_stratum(synthetic, tmp_path):
    """With strata supplied, each window carries its dominant callability tier
    (what the scatter colors by), and the strata + scatter outputs are written."""
    easy = tmp_path / "easy.bed"
    easy.write_text("chr1\t0\t20\n")            # 20 bp of the 30 bp chr1
    diff = tmp_path / "difficult.bed"
    diff.write_text("chr1\t20\t30\n")           # the other 10 bp
    cfg = synthetic.model_copy(update={"strata": {"easy": str(easy), "difficult": str(diff)}})

    rep = QCReport(cfg)
    rep.run(bedgraph_dir=tmp_path / "perbase", jobs=1)

    assert all("stratum" in w for w in rep.win_rows)
    chr1 = [w for w in rep.win_rows if w["chrom"] == "chr1"]
    assert chr1 and chr1[0]["stratum"] == "easy"   # 20 bp easy dominates 10 bp difficult

    written = rep.write_outputs(tmp_path)
    assert written["scatter"].exists()
    assert written["strata"].exists()


def _write_fake_strata(directory):
    """The three SMaHT-named BEDs covering the synthetic contigs."""
    import gzip
    directory.mkdir(parents=True, exist_ok=True)
    beds = {
        "SMaHT_easy_hg38.bed.gz": "chr1\t0\t20\n",
        "SMaHT_difficult_hg38.bed.gz": "chr1\t20\t30\n",
        "SMaHT_extreme_hg38.bed.gz": "chr2\t0\t5\n",
    }
    for name, content in beds.items():
        with gzip.open(directory / name, "wt") as fh:
            fh.write(content)


def test_cli_strata_flag(synthetic, tmp_path):
    """`--strata` is a boolean flag: it picks up the fixed SMaHT tier set from
    --strata-dir (no label=path needed) and writes the per-tier table."""
    sdir = tmp_path / "strata"
    _write_fake_strata(sdir)
    outdir = tmp_path / "out"
    res = CliRunner().invoke(main, ["coverage", "--cram", str(synthetic.cram),
                                    "--reference", str(synthetic.reference), "--full",
                                    "--outdir", str(outdir),
                                    "--strata", "--strata-dir", str(sdir)])
    assert res.exit_code == 0, res.output
    assert (outdir / "coverage.strata.tsv").exists()


def test_cli_strata_flag_errors_without_beds(synthetic, tmp_path):
    """--strata with no BEDs present points the user at `fetch strata`."""
    res = CliRunner().invoke(main, ["coverage", "--cram", str(synthetic.cram),
                                    "--reference", str(synthetic.reference), "--full",
                                    "--outdir", str(tmp_path / "out"),
                                    "--strata", "--strata-dir", str(tmp_path / "nope")])
    assert res.exit_code != 0
    assert "fetch strata" in res.output
