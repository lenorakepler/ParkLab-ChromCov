# Handling file output location in a CLI tool like `chromcov`

Notes on where a coverage tool should read from and write to, especially when run outside the project root or inside a workflow engine. The current `__main__` blocks in `chromcov/*.py` use `Path().resolve() / "data"`, which silently assumes the working directory is the project root; run from anywhere else and it breaks. That is the fragility this document addresses.

## The core principle: separate "where the user is" from "where the package lives"

Almost every output-location bug comes from conflating these two. There are three distinct location sources, and each has exactly one correct resolution method.

| Kind of path | Anchor it to | How |
|---|---|---|
| **User inputs** (CRAM, reference) | the caller's CWD, or absolute | `Path(arg).expanduser().resolve()` — never relative to the package |
| **User outputs** (coverage table) | CWD by default, or an explicit `-o` / prefix | see below |
| **Bundled package data** (test fixtures) | the *package*, not CWD | `importlib.resources` or `Path(__file__).parent` |

The current `data/` lookup does the third thing (finding bundled files) with the first thing's method (CWD-relative), which is the whole problem. If `data/` holds dev fixtures, resolve it via `__file__`; if it holds user data, it should be a required CLI argument.

## Output-location conventions (pick one, don't invent)

Two dominant patterns in this ecosystem:

1. **Output-prefix model** — for tools that emit *multiple related files*. mosdepth is canonical: pass a prefix and it writes `PREFIX.mosdepth.summary.txt`, `PREFIX.per-base.bed.gz`, etc. samtools, plink, and GATK use variants. If `chromcov` grows to emit a summary plus per-base plus breadth table, this is the right convention.
2. **Single `-o FILE`, default stdout** — for a single result. `samtools coverage` and `samtools depth` default to stdout so they are pipeable. Primary result to stdout, logs and progress to stderr, and let the user redirect. For a one-table coverage tool this is the most composable choice.

Rules that hold across both:

- Default output goes to CWD (or stdout), never the install directory. The user expects files where they are standing.
- Prefixes and paths are resolved relative to CWD, and the output directory is created with `mkdir(parents=True, exist_ok=True)` if a nested prefix is given.
- Do not overwrite silently when writing named files — a `--force` flag or a refuse-and-warn is typical.
- Prefer atomic writes (write to a temp file, then `os.replace`) so a killed run does not leave a truncated file that looks complete.
- Temp and intermediate files go through `tempfile` respecting `$TMPDIR`, or a `--tmpdir` flag — never dumped in CWD or the package directory.

## Why this matters specifically for this role: workflow engines

This tool is destined for CWL (per the job description), and that hard-constrains output handling. Under CWL, WDL, Nextflow, and Snakemake, the engine creates a fresh working directory per task, stages inputs into it, runs the tool with that directory as CWD, and collects outputs by globbing filenames in that directory. CWL's `outputBinding: { glob: "*.summary.txt" }` literally picks up files by pattern from the task's CWD.

The implication is strong and specific: a CWL-friendly tool must write to CWD (or an engine-provided path) with predictable, declarable filenames, and must never write to an absolute path or a package-relative path. If the tool hardcoded an output location, the engine could not capture the result at all. So "default to CWD, name outputs predictably (ideally via a prefix)" is not just tidy — it is the thing that makes the tool composable into their pipelines, and worth a sentence in the writeup.

## Concrete shape for `chromcov`

Roughly what the real entry point (as opposed to the `__main__` smoke tests) would do:

```python
# inputs: explicit, resolved absolutely, no package/CWD-relative assumptions
cram = Path(args.cram).expanduser().resolve()
ref  = Path(args.reference).expanduser().resolve()

# output: default to stdout; -o writes to a CWD-relative path
if args.output:
    out = Path(args.output).expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    fh = out.open("w")
else:
    fh = sys.stdout            # pipeable default

# bundled test fixtures (if any): located via the package, not CWD
from importlib.resources import files
fixture = files("chromcov") / "data" / "tiny.cram"
```

## One-line takeaway

Inputs and outputs anchor to the user's CWD (explicit, absolute-resolved, stdout by default); only bundled package data anchors to the package via `__file__` / `importlib.resources`; and writing to CWD with predictable names is precisely what lets a CWL engine capture the results. Fixing that in the `__main__` blocks is the difference between a script and a tool.
