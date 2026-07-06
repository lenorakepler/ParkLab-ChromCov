# ChromCov

This repo is a general-purpose tool for calculating the average per-base sequence coverage for each chromosome in a CRAM file.

## Installation

```
  git clone https://github.com/lenorakepler/ParkLab-ChromCov.git
  cd ParkLab-ChromCov
 
  # builds the env from uv.lock + installs the `chromcov` CLI
  uv sync                
```

No `uv`? `pip install -e .` works too.

## Usage

By default, the the coverage calculation excludes unmapped, secondary, and supplementary reads. Since COLO829T is a cancer cell line, does not filter out duplicated reads. Overlapping mate pairs handled by the `pysam.AlignedSegment.get_blocks()`

### Quick, no-frills

At its simplest, one can use program defaults and just point to the CRAM and reference files. The program will output per-chromosome mean coverage to stdout.

```bash
# install the package and CLI
uv sync

chromcov coverage \
  --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

# write it to a file instead
chromcov coverage --cram … --reference … -o out/coverage.tsv
```

### A little more fun (`--full`)

For more context around the average-coverage number, run the program with `--full`

```bash
chromcov coverage --cram … --reference … --full --outdir out/ --jobs 4
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

| path                    | description                                                           |
| ----------------------- | --------------------------------------------------------------------- |
| `chromcov/coverage.py`  | Creates per-chromosome mean coverage table                            |
| `chromcov/calc_cov.py`  | The coverage calculation itself (event-based, O(reads))               |
| `chromcov/config.py`    | Config, including read flags                                          |
| `chromcov/report.py`    | Extended statistics report tables                                     |
| `chromcov/qc_report.py` | `--full` orchestration                                                |
| `chromcov/depth.py`     | per-base reductions: `ChromDepth` · `DepthHistogram` · `ChromStats`   |
| `chromcov/qc_flags.py`  | QC abnormality flags (user-defined values in config)                  |
| `chromcov/strata.py`    | SMaHT easy/difficult/extreme callability tiers + masks                |
| `chromcov/plots.py`     | Bar and scatter coverage plots                                        |
| `chromcov/perbase.py`   | Generate bedgraphs + resume                                           |
| `chromcov/validate.py`  | Ensure CRAM-Reference validity                                        |
| `chromcov/fetch.py`     | download the SMaHT strata BEDs                                        |
| `chromcov/cli.py`       | CLI (`coverage` / `plot` / `fetch`) +`run.json`                       |
| `tests/`                | pytest: reduction/QC math + a synthetic-CRAM end-to-end coverage test |

## Write-Up

The main thing is a custom function to calculate coverage (calc_cov.py:calc_cov). It uses the same algorithm as mosdepth, where instead of recording and summing every read at every position, O(depth * reference positions) it instead records the start and end of each read, ~O(reads). It adds +1 to the first reference position a read covers a -1 to the first position it does *not* cover. A cumulative sum (truly it's like magic) yields the total coverage at each base in the reference genome.

I tried to look at where this metric might be integrated into existing workflows by examining both the Park Lab repo and the SMaHT DAC repo. This confirmed my suspicion that per-chromosome coverage is a broad QC metric that is unlikely to be used as input into any downstream analyses. This statistic can be used to confirm at larger scales to confirm that target depth is adequate and on a per-sample basis to surface gross, large-scale feature estimates like aneuploidy and the sex-chromosome complement that might themselves be verifying QC signals or, alternately, previously-unknown features that might bias results or be worthy of further investigation.

However, I also wanted to provide additional extra context to that number and give the option to get a report of the overall shape of the coverage distribution (MAD, IQR); how coverage relates to reference position, to allow for quick visualization of potential breakpoints; coverage breadth at a given sequencing depth; and comparison between chromosomes and to an average autosomal depth to use as a signal of possible aneuploidy.

Everything is wrapped with some machinery: a config to parse options that can be overridden with the CLI; json output documenting run specifications and code state to ensure reproducibility; validation to ensure that the correct genome reference is being used and that the CRAM is sorted, which is necessary for the coverage calculation; and a function to pull down files from AWS. These may seem superfluous, but reproducibility and ease-of-use seemed important to the "helpful for other bioinformatics analysts and software engineers" piece.

I also noticed the "SMaHT_Regional_Categorization" repo and thought it would be fun to integrate those metrics into some of the plots, to examine to what extent the coverage calculation is being biased by sequencing artifacts, etc.

This can be threaded at the command line, or used in a workflow manager like CWL. It's packaged with UV, which I have found to lead to an extremely smooth installation and run experience that will not mess with any existing python environments.

I used Claude Code to help in writing some peripheral functions like the flag filter mask and code tests, as well as in sketching out the architecture and class structure I wanted over many iterations and helping to debug generally.
