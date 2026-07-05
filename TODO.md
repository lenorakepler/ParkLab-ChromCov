# TODO — path to a finished, sendable repo

Ordered by what a reviewer will notice first. Items marked **(blocking)** should
be done before sending; the rest strengthen the submission.

## Scope & framing (read this first — my honest pushback)

The literal ask is *"average per-base coverage per chromosome."* That core is a
few lines: `aligned_bases / chromosome_length`. Everything else here (two
backends, strata, copy number, QC flags, plots) is genuinely useful and shows
range — **but it must be framed as deliberate extension, not mistaken for the
deliverable.** Two risks to manage:

1. **Don't let the extras bury the core.** The README leads with the plain answer
   on purpose. In the written summary, state the simple result and method first,
   then "beyond the ask, here's what a production QC step needs." Over-engineering
   without that framing reads as poor scope judgment — the opposite of the signal
   you want.
2. **Decide how much to ship.** Given it's an interview, I'd lean toward: a *crisp,
   installable core package* + the analysis suite clearly labeled as extensions +
   a strong written summary. Shipping the current `dev/` sketches as-is is the
   thing I'd most push back on (see next section).

## Blocking — consolidate sketches into the package

Right now all the polished code is runnable **sketches** under
`dev/dispatch-sketch/`, while the actual `chromcov/` package is empty/broken. A
reviewer will `pip install` and run the package, not read `dev/`.

- [ ] **(blocking)** Move the sketch modules into `chromcov/` proper
      (`analysis`, `strata`, `qc`, `plots`, `output`, `validate`, backends).
- [ ] **(blocking)** Delete the `get_cov.py` / `read_filter.py` duplication — keep
      the `ReadFilter` version; drop `get_cov.py`.
- [ ] **(blocking)** Fix or delete `chromcov/stats.py` (references an unimported
      `np`, returns undefined names — superseded by `analysis.py`).
- [ ] **(blocking)** Add a CLI entry point: `[project.scripts]
      chromcov = "chromcov.cli:main"` with subcommands `coverage` / `analyze`.
      The CWL tool and Dockerfile already assume this contract.
- [ ] **(blocking)** Wire `config.example.yml` to a loader (`CoverageConfig.from_yaml`)
      so a run is driven by one file (needs `pyyaml`).

## Data acquisition (requested)

- [ ] **Add a `chromcov fetch-data` command / script** that downloads the three
      inputs from the URLs in the assignment email into `data/`:
      - `COLO829T_TEST.cram`
      - `COLO829T_TEST.cram.crai`
      - `GCA_000001405.15_GRCh38_no_alt_analysis_set.fa`
      (base: `https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/`)
      Verify sizes/checksums after download and skip files already present. This
      makes the repo runnable from a clean clone and documents provenance of inputs.
- [ ] **Support reading the CRAM directly from a URL** where feasible.
      Notes for whoever implements it: htslib/pysam *can* open a remote CRAM and
      fetch by region over HTTP range requests using a (local or remote) `.crai` —
      `pysam.AlignmentFile("https://…/…cram", index_filename=…)`. The catch is the
      **reference**: CRAM decode needs the reference available (locally or via
      `REF_CACHE`/`REF_PATH`), and random reference access over the network is slow,
      so the reference realistically still has to be downloaded once. So: stream
      the CRAM from S3, but stage the reference locally (or a `REF_CACHE`). Make
      remote-vs-local an input mode, not the only path.
- [ ] Generate + commit a `samtools dict` for the reference (enables the fast,
      no-rehash M5 verification in `validate.py`).

## Tests & CI

- [ ] **(blocking)** Unit tests against a **tiny synthetic CRAM** with hand-computed
      known depth: coverage math (incl. deletions/skips contributing 0), aggregate
      vs per-base agreement, `stats_from_hist` (mean/median/MAD/breadth), windowing,
      `stratum_mask`, and `qc` flags.
- [ ] Run the native↔mosdepth cross-check (`test_backends.py`) in CI on Linux,
      where mosdepth installs cleanly.
- [ ] GitHub Actions workflow: lint (ruff), tests on the fixture, and a Docker build.
      The JD calls out CI/CD explicitly.

## Correctness / robustness

- [ ] Make `LOW_CALLABLE` (and other depth thresholds) **baseline-relative** rather
      than the fixed 20× that over-flags a 15× sample.
- [ ] Add a `--no-overlap` option to avoid double-counting overlapping mate pairs
      (document that the default matches mosdepth, which also double-counts).
- [ ] Run the **mosdepth backend end-to-end on Linux** to confirm parity (never
      completed here due to local Docker/macOS-FUSE trouble).
- [ ] Parallelize the native backend by chromosome (`multiprocessing`) if the ~7 min
      single-threaded genome pass needs to be faster; or lean on mosdepth for speed.

## Reproducibility / packaging

- [ ] Fold `preflight`'s report into the provenance sidecar (it already records
      reference verification — merge the sorted/index results too).
- [ ] Add `mosdepth` to the reproducibility Dockerfile so `docker run chromcov`
      exposes both backends in one image; publish the image for the CWL tool.
- [ ] Pin plotting deps (`matplotlib`) in `pyproject.toml` rather than
      `uv run --with matplotlib`.
- [ ] BigWig output option for per-base (pyBigWig / bedGraphToBigWig) — see
      DEVELOPMENT §7 on why flat per-base bedgraph is a poor default.

## Nice-to-have

- [ ] Per-chromosome × stratum breakdown (currently strata are pooled genome-wide).
- [ ] HTML/PDF report bundling the table, plots, and flags for the "analyst" audience.
- [ ] GC-bias normalization for the copy-number estimate.
