# Dispatch sketch

A sketch of the two-backend coverage dispatcher: one shared `CoverageConfig` in,
one normalized `list[ChromCoverage]` out, backend (hand-rolled pysam vs.
`mosdepth`) chosen by config. Illustrative — not yet wired into `chromcov/`.

| File | Role |
|------|------|
| `config.py` | `CoverageConfig` — the shared knobs; pins an explicit flag mask so both backends agree. |
| `result.py` | `ChromCoverage` — one row per chromosome (`mean = bases / length`), plus a TSV writer. |
| `mosdepth.py` | subprocess wrapper + `.mosdepth.summary.txt` parser. |
| `dispatch.py` | `run_coverage(config)` — picks the backend. |
| `native.py` | adapter over `chromcov.read_filter.calc_cov`. |
| `output.py` | per-run archival (hash-keyed dir + sidecar) and cross-run `collate`. |
| `validate.py` | input preflight: sorted / indexed / reference-M5 match. |
| `run_example.py` | end-to-end smoke run against `data/`, with `--write` / `--collate`. |
| `test_backends.py` | native-vs-mosdepth cross-check (skips without mosdepth/data). |
| **analysis suite** | |
| `analysis.py` | histogram → stats (mean/median/sd/CV/MAD/breadth), windows, RLE, copy number. |
| `strata.py` | SMaHT easy/difficult/extreme callability loader + finite-diff mask. |
| `qc.py` | aneuploidy-aware abnormality flags + focal per-window flags. |
| `plots.py` | per-chrom bar + windowed-CN scatter (easy-maskable). |
| `analyze.py` | driver → stats/windows/strata TSVs, plots, optional per-base bedgraph. |
| `config.example.yml` | full annotated run-config sketch. |

## The native backend adapter (not written here)

`native.run` is a thin shim over the existing `chromcov.read_filter.calc_cov`.
It builds a `ReadFilter` from the shared config and loops chromosomes:

```python
# native.py (sketch)
import pysam
from chromcov.read_filter import ReadFilter, calc_cov
from config import CoverageConfig
from result import ChromCoverage

def run(config: CoverageConfig) -> list[ChromCoverage]:
    rf = ReadFilter(
        include_flags=config.include_flags,
        exclude_flags=config.exclude_flags,
        exclude_all_flags=config.exclude_all_flags,
        min_mapping_quality=config.min_mapping_quality,
    )
    cram = pysam.AlignmentFile(
        str(config.cram), "rc",
        reference_filename=str(config.reference),
        index_filename=str(config.index),
    )
    rows = []
    for chrom in cram.references:
        base_depth, total_depth, _ = calc_cov(cram, chrom, rf, per_base=config.per_base)
        rows.append(ChromCoverage(
            chrom=chrom,
            length=cram.get_reference_length(chrom),
            bases=int(total_depth),
            backend="native",
            base_depth=base_depth,
        ))
    return rows
```

(When merging `get_cov.py` into `read_filter.py`, drop the stray `breakpoint()`
in `calc_cov`.)

## Why the numbers should match

Both count aligned bases over M/=/X CIGAR ops (mosdepth normal mode == your
`get_blocks()`), neither corrects for overlapping mates, and `mean` is
`bases / length` in both. With `exclude_flags` pinned to the same mask in
`config.py`, the two backends become a mutual check:

```python
n = run_coverage(CoverageConfig(cram, ref, backend="native"))
m = run_coverage(CoverageConfig(cram, ref, backend="mosdepth"))
# assert per-chrom means agree to rounding
```

`test_backends.py` encodes exactly that (on a single small contig via mosdepth
`-c`), and skips when mosdepth or `data/` is absent:

```bash
uv run --with pytest pytest dev/dispatch-sketch/test_backends.py -v
```

## Output handling: archive per run, collate to compare

Two jobs, two shapes (`output.py`):

- **Archival** — `write_run` puts each run in `runs/<name>/`. Inside: a
  self-describing TSV (`#` comment header with tool version, commit, params,
  input paths) sorted karyotypically, plus the `provenance.json` sidecar from
  `reproducibility-sketch/provenance.py`. Two naming styles (`--run-name`):
  - `slug` (default) — readable `backend[-deviations]-<hash>`, e.g.
    `native-q20-a8654603`. Shows only params that differ from default
    (samtools letters `q`/`f`/`F`/`G`, masks in hex), so `ls runs/native-q20-*`
    just works.
  - `hash` — the bare digest, e.g. `a8654603`.

  Both end in the same deterministic hash of `(params + input identity)`, so
  either way: same config → same dir → idempotent (re-running is a detected
  no-op); different config → different dir → never clobbers. The hash suffix is
  what makes the friendly-but-lossy slug still unique (it disambiguates params
  and inputs the slug hides).
- **Comparison** — `collate` walks `runs/*/`, reads each sidecar for its params,
  and stacks everything into one long-format table; `pivot_mean` reshapes to wide
  (chrom × run). So comparing backends or filters is a group-by, not a manual diff.

```bash
uv run python dev/dispatch-sketch/run_example.py --chroms chr21 --min-mapq 0  --write
uv run python dev/dispatch-sketch/run_example.py --chroms chr21 --min-mapq 20 --write
uv run python dev/dispatch-sketch/run_example.py --collate
# chrom   native-6853f9b2   native-q20-a8654603
# chr21   16.65             12.65     <- MAPQ 20 drops low-quality alignments
```

(`runs/` is gitignored. The header shows `vNone` until the package is installed
so `importlib.metadata` can resolve its version — a real install via
`[project.scripts]` fixes that; the git commit + dirty flag work regardless.)

## Calling mosdepth

`mosdepth.py` calls the binary on PATH — not `docker run`. The binary is pinned
by the environment (add `mosdepth` to the project Dockerfile via the biocontainer
or bioconda), which is what the CWL runner already containerizes. Locally:

```bash
mamba install -c bioconda mosdepth   # or grab the static release binary
```
