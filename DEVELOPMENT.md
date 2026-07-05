# Development notes

Design rationale and the details that aren't obvious from the code. Aimed at the
next engineer (or the reviewer) who needs to understand *why*, not just *what*.

## 1. The coverage algorithm

Naïve per-base coverage increments a counter at every position each read covers —
O(coverage × genome). Instead we accumulate **events**: for each aligned block
`[start, end)`, do `depth[start] += 1` and `depth[end] -= 1`. A single `cumsum`
at the end turns that finite-difference array into per-base depth. Cost is
O(number of alignment blocks) ≈ O(reads), independent of depth.

Two consequences worth stating explicitly:

- **The arithmetic is order-independent.** `+1`/`-1` then `cumsum` is commutative,
  so the tally does not depend on read order. Coordinate-sorting is therefore *not*
  required by the coverage math — it's required by the **indexed region fetch**
  (`fetch(chrom)` needs a `.crai`, which only indexes a sorted CRAM). See
  `validate.py`; this is the kind of thing that's easy to state wrongly.
- **Block semantics matter.** We use `read.get_blocks()`, which yields the gapless
  aligned blocks defined by CIGAR `M`/`=`/`X` (match/mismatch), excluding `I`/`D`/
  `N`/`S`/`H`. So deletions and reference skips correctly contribute **no** depth —
  the same convention mosdepth uses in normal mode, which is why the two backends
  agree to rounding.

For the aggregate-only path (`per_base=False`) we skip the array entirely and just
sum `end - start` — enough for the headline mean, in less memory.

## 2. Memory strategy

The per-base vector for chr1 is ~1 GB (int32 × 249 Mb). We never hold the whole
genome: `CoverageAnalysis` (`pipeline.py`) processes **one chromosome at a time**,
wrapping each per-base vector in a `ChromDepth` (`analysis.py`), immediately
reducing it to compact intermediates, then `del`s it. Peak memory ≈ one chromosome.

The two intermediates are chosen because each is a *sufficient statistic* for a
whole family of outputs:

- **Depth histogram** (`ChromDepth.histogram`) — counts per depth value. Yields
  mean, median, variance, CV, MAD, any quantile, and breadth-at-depth, all in
  O(max_depth) with no sort. Genome-wide stats = summed per-chrom histograms.
- **Windowed mean track** (`ChromDepth.windowed_means`, `np.add.reduceat`) — ~1000×
  smaller than per-base; feeds every plot and the per-window copy number.

Extreme pileup depths (chrM ~20k, satellite decoys ~2k) are clipped into a top
histogram bin (`DEFAULT_HIST_CAP`) so per-chrom histograms share a length and sum
trivially. This is exact for the primary assembly; only the extreme tail of
artifact contigs is approximated.

## 3. Robust statistics

Multi-mapping pileups give coverage a heavy right tail, which wrecks moment-based
stats: chr21 has `sd = 71`, `cv = 4.25` — meaningless. So we report **median**,
**scaled MAD** (`1.4826 × median(|x − median|)`, a consistent sd estimator under
normality), and **robust CV = MAD/median** alongside the classical ones. The
`UNEVEN` QC flag keys off robust CV, so a repeat pileup doesn't trigger a false
"non-uniform coverage" alarm. MAD is computed as a weighted median of folded
deviations, straight from the histogram.

## 4. One calculator, mosdepth as an optional cross-check

There is a single in-tool calculator: the native `calc_cov` (`dispatch.run_coverage`).
mosdepth is **not** a second backend baked into the dispatch path — it's an
optional add-on, `scripts/mosdepth_coverage.py`, that runs the `mosdepth` binary
and converts its output into chromcov's `coverage.tsv` format. So validation is
an explicit, opt-in `diff` (`tests/test_mosdepth_compare.py`) rather than a
registry/ABC/subprocess seam the core has to carry. It also keeps the deliverable
path dependency-light (no `mosdepth` needed to run chromcov).

Getting parity right still required reconciling **default flag masks**: the
hand-rolled default excluded unmapped|secondary|dup|**supplementary** (3332);
mosdepth's default excludes unmapped|secondary|**qcfail**|dup (1796). `config.py`
pins the explicit union (3844) as `DEFAULT_EXCLUDE`, and the add-on defaults to
that same mask, so the two agree out of the box. mosdepth has no `-G`
(exclude-all-flags) equivalent, so the add-on simply doesn't offer that knob.

The add-on shells out to `mosdepth` on PATH, never via `docker run` from Python:
containerization belongs at the environment boundary (the CWL/Docker image that
wraps the whole tool), not nested inside it — nested containers mean
docker-in-docker on AWS/HPC and host↔container path translation.

## 5. Callability stratification

Coverage without callability context is partly an artifact: raw per-chromosome
means are inflated by repeat pileups and deflated by unmappable regions. We
stratify by Park Lab's **SMaHT_Regional_Categorization** tiers — `easy` (1000G
strict mask), `difficult` (PanMask pm151 minus 1000G), `extreme` (outside both) —
loaded from their committed BEDs into a `Strata` (`strata.py`). `Strata.mask`
turns each stratum into a boolean position mask using the *same* finite-difference
trick as the coverage calc, then `ChromDepth.masked(mask).histogram()` reduces the
per-base vector into a per-tier histogram.

Payoff: `easy` mean/median is the **variant-callable** coverage worth reporting,
`easy.breadth_20x` is the somatic-sensitivity metric, and the CN baseline uses
only easy-autosomal positions so repeat pileups don't inflate the CN=2 reference.

## 6. Approximate copy number

`CN ≈ ploidy × depth / baseline`, where `baseline` is the **median** (robust)
depth of the callability-masked autosomes. Median (not mean) so gains/losses and
pileups don't drag the diploid reference. Reported per chromosome and per window;
the windowed scatter (callable windows only) shows intrachromosomal breakpoints.

Deliberately labeled **super approximate**: it ignores tumor purity, ploidy
normalization, and GC/mappability bias. Real callers (CNVkit, ichorCNA) handle
those; this is a QC-grade readout, not a CN call.

## 7. Per-base output — what we learned

The finite-difference change-points *are* a run-length encoding, so the per-base
depth tracks (`PerBaseStore`, via `ChromDepth.rle_intervals`) are RLE BEDGRAPH
rather than one row per base. **But** on dense WGS this only compresses ~6× (chr21:
7.07M intervals for 46.7M bases; chrM changes almost every base), because depth
flips at every read boundary. So the tracks are the reusable base layer but are
bulky; `--fast` skips them entirely (mean only), and **BigWig** is the better
format if a genome-browser track is truly needed (its zoom levels compress far
better than flat bedgraph).

## 8. Preflight & provenance

`validate.py::preflight` fails fast, cheapest-first: inputs exist → `@HD
SO:coordinate` → `.crai` present → **reference M5 matches the CRAM**. The M5 check
is the important one: CRAM stores bases as differences against the reference, so
the wrong reference can silently reconstruct wrong bases. We compare the CRAM's
`@SQ M5` tags to the reference's per-sequence MD5, preferring a `.dict` sidecar
(no re-hashing) over reading the 3 GB `.fa`. `provenance.py` records code commit +
params + input identity + this reference verification in a `*.provenance.json`
sidecar, so every output traces back to an exact run.

A run is assembled in one place — `RunConfig.load` (config.py): the `--config`
file is the base, the CLI contributes only overrides. The output sidecars then
embed the **resolved** `RunConfig` (plus a pointer to the source config file), so
a per-base track (`coverage.json`) or an analysis run (`run.json`) is
self-describing and re-runnable from the sidecar alone.

## 9. Testing approach

`tests/test_backends.py` cross-validates native vs mosdepth on one contig (skips
when mosdepth/data are absent). `tests/test_analysis.py` unit-tests the reduction
math directly on hand-built depth vectors — `ChromDepth` histogram/stats
(incl. the lower-quantile convention)/windowing/RLE, `Strata.mask` (overlap +
clip), the aneuploidy-aware `qc` flags, and copy number — needing no CRAM. The
remaining gap (see `TODO.md`) is a tiny synthetic CRAM to exercise `calc_cov`
end-to-end (deletions/skips contributing 0 depth) against known values.

## Known correctness caveats

- **Overlapping mate pairs are double-counted.** `get_blocks` doesn't dedupe
  overlapping read1/read2; mosdepth's default doesn't either, so the backends
  agree — but both overstate depth in overlaps. A `--no-overlap` option is a real
  follow-up.
- **`LOW_CALLABLE` threshold is absolute (20×).** It should be baseline-relative
  (e.g. breadth at ~½ the diploid baseline); on a 15× sample the fixed 20× check
  over-flags.
- **Histogram cap** approximates the extreme tail of pileup contigs (not the
  primary assembly).
- **mosdepth backend is unrun end-to-end** here (local Docker/macOS-FUSE trouble);
  parity is by construction and needs a Linux confirmation run.
