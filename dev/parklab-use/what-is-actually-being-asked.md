# What the take-home is actually asking for

Companion to [`coverage-downstream-research.md`](./coverage-downstream-research.md). This one steps back from *how coverage is computed downstream* to *what this assignment is really testing*, and why its vagueness is deliberate.

## The prompt, and its two halves

> "Develop a program that calculates the average per-base sequence coverage for each chromosome in the CRAM file we are sharing. Provide your code and report the answer in ways you think are optimal/helpful for other bioinformatics analysts and software engineers in a team."

The first sentence describes a computation anyone on the team could run in one line: `samtools coverage COLO829T_TEST.cram`. If they wanted the number, they'd already have it. So the operative brief is the second sentence — reporting it "in ways optimal/helpful for other bioinformatics analysts and software engineers." This is a judgment, communication, and reproducibility exercise wearing a trivial-computation costume, and it maps almost line-for-line onto the job description (reproducible code, CWL, Docker, AWS, testing, documentation).

## Why it's underspecified on purpose

"Average per-base coverage per chromosome" sounds precise but leaves several decisions open, and every one of them changes the number. The test is whether the candidate surfaces these decisions and makes defensible, documented choices rather than picking silently. The decision points:

1. **What counts as "depth"?** Per-base pileup depth versus alignment-block coverage, and — critically — whether overlapping mate pairs are counted once or twice. This is exactly the mosdepth default-versus-`--fast-mode` distinction: default mosdepth corrects mate overlaps and is CIGAR-aware, while `--fast-mode` disables both. Notably the SMaHT DAC itself runs `mosdepth -n --fast-mode`, i.e. it chose the faster, double-counting version. There is no single right answer here, only a *stated* one. (For reference, `chromcov` in this repo is CIGAR-aware via `read.get_blocks()` like default mosdepth, but does not correct mate overlaps like fast-mode mosdepth — a hybrid, which is fine as long as it's documented.)

2. **Which reads count?** Duplicates, secondary and supplementary alignments, unmapped reads, a MAPQ floor, a base-quality floor. Real pipelines disagree: the DAC's Picard `CollectWgsMetrics` uses MQ/BQ ≥ 20, Park Lab's SCAN2 uses MQ 60, and `chromcov`'s default drops unmapped/secondary/supplementary/duplicate with MQ 0. All are defensible; all yield different numbers.

3. **What is the denominator?** Averaging over the full chromosome length (zero-depth bases included) versus over only covered bases. For a cancer genome with hemizygous deletions these diverge substantially. "Average per-base coverage" implies full length, but the choice should be explicit.

4. **What is a "chromosome"?** The `GCA_000001405.15_GRCh38_no_alt_analysis_set` reference has roughly 195 contigs. chrM reads absurdly high, decoy/EBV/unplaced contigs near zero. Report all of them, only the primary 1–22/X/Y, or grouped? Any of these is fine if stated.

5. **The sample context.** COLO829T is an aneuploid cancer cell line with many copy-number variations, so per-chromosome means *should* deviate from the genome-wide average — that is biological signal, not an error. Noticing this signals understanding of the data, not just the file format.

## What "done" looks like, in tiers

**Table stakes (necessary, not sufficient):** open the CRAM with the provided reference, emit per-chromosome mean depth, and land within tolerance of `samtools coverage` / mosdepth.

**What actually earns the interview:**

- **Explicit, documented methodology** — the decision points above, written down with rationale.
- **Validation** against a standard tool, showing the number is close to mosdepth / samtools and explaining any delta (e.g. mate-overlap correction).
- **Reproducibility** — a pinned environment plus Docker, and ideally a CWL wrapper. The job description asks for exactly this, so handing them a CWL tool is speaking their language.
- **Tests** — even a couple against a tiny fixture.
- **Dual output** — a human-readable report *and* machine-readable TSV/JSON. This is the "analysts *and* engineers" clause taken literally: analysts want the interpreted table plus caveats, engineers want parseable output plus a clean CLI.
- **A short written summary** — choices, caveats, the CNV interpretation, and a note on scale (this is a subsampled test CRAM; real WGS is hundreds of GB, so how would it parallelize / run on AWS).

## The one-sentence reframe

The tool needs to compute a number anyone could get in one line, so the deliverable is not the number — it is the evidence that the candidate makes careful methodological choices and ships reproducible, well-communicated software. The vagueness is not a gap to fill silently; it is the surface being evaluated, so every hidden decision should be made visible.
