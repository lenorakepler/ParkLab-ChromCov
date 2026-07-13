# ChromCov

This repo is a general-purpose tool for calculating the average per-base sequence coverage for each chromosome in a CRAM file.

## Installation

```
  git clone https://github.com/lenorakepler/ParkLab-ChromCov.git
  cd ParkLab-ChromCov
 
  # builds the env from uv.lock + installs the `chromcov` CLI
  uv sync

  # activate the env so the bare `chromcov` command is on your PATH
  source .venv/bin/activate
```

No `uv`? `pip install -e .` (into an activated venv) works too.

Prefer not to activate? Prefix any command below with `uv run`, e.g. `uv run chromcov coverage --config config.yaml`. To put `chromcov` on your PATH globally instead, `uv tool install --editable .` (note: this resolves fresh and does not use `uv.lock`).

## Getting the data

The tool defaults to the Park Lab COLO829T test files under `data/`. If your clone doesn't already have them, download them (the CRAM + reference are large -- tens of GB):

```bash
chromcov fetch inputs   # -> data/COLO829T_TEST.cram (+ .crai) and the GRCh38 reference
```

Once fetched, a bare `chromcov coverage` runs on them with no flags. Point `--cram` / `--reference` (or the config's `inputs:` block) elsewhere to run on your own data.

## Usage

By default, the the coverage calculation excludes unmapped, secondary, and supplementary reads. Since COLO829T is a cancer cell line, does not filter out duplicated reads. Overlapping mate pairs handled by the `pysam.AlignedSegment.get_blocks()`

### Recommended options

It is highly recommended to generate and use a config file. With this command, one can generate a config file that exposes all program defaults so that you can make sure it works as you are expecting and re-configure anything you need:

```bash
chromcov gen-config [-o config.yaml]
```

### Quick, no-frills

At its simplest, one can use program defaults, which point at the bundled COLO829T test files under `data/`. The program will output per-chromosome mean coverage to stdout.

```bash
# runs on the default data/ inputs -- no flags needed
chromcov coverage

# point at your own CRAM + reference
chromcov coverage \
  --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# or drive everything from a config file
chromcov coverage --config config.yaml
```

### More-context QC (`--full`) -- a little more fun!

For more context around the average-coverage number, run the program with `--full`. If analyzing the whole genome, chromosomes can be analyzed in parallel with `--jobs`

```bash
chromcov coverage --config config.yaml --full --outdir out/ --jobs 4
```

This writes, under `out/`:

| file                                | description                                                                             |
| ----------------------------------- | --------------------------------------------------------------------------------------- |
| `perbase/chrN.per-base.bedgraph.gz` | per-base depth, one per chromosome                                                      |
| `coverage.tsv`                      | extended coverage statistics table: mean, median, IQR, MAD, breadth-at-deapth, QC flags |
| `coverage.windows.bed`              | windowed mean depth + copy number + focal flag along the genome                         |
| `coverage.bar.png`                  | mean depth per chromosome                                                               |
| `coverage.scatter.png`              | windowed copy number, colored by callability tier                                       |
| `run.json`                          | the resolved config + code state, so the run reproduces from the output dir alone       |

### Other options

#### Specify a config file

`--config run.yaml`

#### Filter reads by mapping quality

`--min-mapq N`

#### Run only a subset of chromosomes:

`--chroms chr21,chrM` (or `21,M` depending on which reference is used)

#### Stratify sites based on how difficult mapping/variant calling is

Download BED files from https://github.com/parklab/SMaHT_Regional_Categorization

1. `chromcov fetch strata`
2. Run `--full` with `--strata`

#### Set window size for finding + visualizing coverage across chromosomes

Run `--full` with `--window N`

## Code files

| path                       | description                                                             |
| -------------------------- | ---------------------------------------------------------------------- |
| `chromcov/kernel.py`       | The coverage calculation itself — `calc_cov`, event-based, O(reads)    |
| `chromcov/filtering.py`    | SAM read-flag vocabulary + per-read `ReadFilter` (`-Q`/`-f`/`-F`/`-G`) |
| `chromcov/reduce.py`       | Per-base reductions: `ChromDepth` · `DepthHistogram` · `ChromStats`    |
| `chromcov/pipeline.py`     | Single orchestrator: source → reduce → accumulate (mean & `--full`)    |
| `chromcov/result.py`       | `RunResult` — the accumulator a run folds into                         |
| `chromcov/policy.py`       | Coverage-QC policy: diploid baseline, copy number, abnormality flags   |
| `chromcov/categories.py`   | SMaHT easy/difficult/extreme callability tiers + masks                 |
| `chromcov/preflight.py`    | Input validation: sorted CRAM, index present, reference matches        |
| `chromcov/io/alignment.py` | CRAM access: context-managed reader, contig listing + lengths          |
| `chromcov/io/codec.py`     | Reference verification (CRAM `@SQ` M5 vs reference MD5)                 |
| `chromcov/io/track.py`     | Per-base depth track I/O (RLE `bedgraph.gz`); the resume boundary       |
| `chromcov/io/fetch.py`     | Download the COLO829T inputs + SMaHT strata BEDs                        |
| `chromcov/config/schema.py`| Pydantic run-config — the sole source of field truth                   |
| `chromcov/config/template.py`| Generate the editable config YAML from the live model (`gen-config`) |
| `chromcov/present/frames.py`| Assemble + write coverage tables (polars); `--full` output orchestration |
| `chromcov/present/plots.py`| Bar and scatter coverage plots (headless Agg)                          |
| `chromcov/present/sidecar.py`| `run.json` provenance sidecar (resolved config + code state)         |
| `chromcov/cli.py`          | CLI (`coverage` / `plot` / `fetch` / `gen-config`)                     |
| `tests/`                   | pytest: reduction/QC math + synthetic-CRAM end-to-end + mosdepth compare |

## Write-Up

The main thing is a custom function to calculate coverage (kernel.py calc_cov). It uses the same algorithm as mosdepth, where instead of recording and summing every read at every position, O(depth \* reference positions) it instead records the start and end of each read, ~O(reads). It adds +1 to the first reference position a read covers a -1 to the first position it does *not* cover. A cumulative sum (truly it's like magic) yields the total coverage at each base in the reference genome.

I tried to look at where this metric might be integrated into existing workflows by examining both the Park Lab repo and the SMaHT DAC repo. This confirmed my suspicion that per-chromosome coverage is a broad QC metric that is unlikely to be used as input into any downstream analyses. This statistic can be used to confirm at larger scales to confirm that target depth is adequate and on a per-sample basis to surface gross, large-scale feature estimates like aneuploidy and the sex-chromosome complement that might themselves be verifying QC signals or, alternately, previously-unknown features that might bias results or be worthy of further investigation.

However, I also wanted to provide additional extra context to that number and give the option to get a report of the overall shape of the coverage distribution (MAD, IQR); how coverage relates to reference position, to allow for quick visualization of potential breakpoints; coverage breadth at a given sequencing depth; and comparison between chromosomes and to an average autosomal depth to use as a signal of possible aneuploidy.

Everything is wrapped with some machinery: a config to parse options that can be overridden with the CLI; json output documenting run specifications and code state to ensure reproducibility; validation to ensure that the correct genome reference is being used and that the CRAM is sorted, which is necessary for the coverage calculation; and a function to pull down files from AWS. These may seem superfluous, but reproducibility and ease-of-use seemed important to the "helpful for other bioinformatics analysts and software engineers" piece.

I also noticed the "SMaHT_Regional_Categorization" repo and thought it would be fun to integrate those metrics into some of the plots, to examine to what extent the coverage calculation is being biased by sequencing artifacts, etc.

This can be threaded at the command line, or used in a workflow manager like CWL. It's packaged with UV, which I have found to lead to an extremely smooth installation and run experience that will not mess with any existing python environments.

I used Claude Code to help in writing some peripheral functions like the flag filter mask and code tests, as well as in sketching out the architecture and class structure I wanted over many iterations and helping to debug generally.
