# Reproducibility sketch

A minimal, deliberately *proportionate* demonstration of reproducible-code
practices around the chromcov coverage step — the kind the Park Lab / SMaHT DAC
role calls out (CWL, Docker, CI/CD, provenance). These files are illustrative
sketches, not yet wired into the build.

## What's here

| File | Role |
|------|------|
| `provenance.py` | Runnable module: stamps every coverage output with a `<output>.provenance.json` sidecar (code version + params + input identity + CRAM/reference verification). |
| `coverage.cwl` | CWL `CommandLineTool` wrapping the coverage step, with `DockerRequirement`. |
| `inputs.example.yml` | Example CWL job file pointing at `data/`. |
| `Dockerfile` | Pinned two-stage image (installs from the committed `uv.lock`). |

## The three things any output must be traceable to

1. **Code version** — exact git commit + a *dirty* flag (`provenance.py::git_provenance`).
2. **Parameters** — captured verbatim into the sidecar (`params=...`).
3. **Environment** — pinned by `uv.lock` and baked into the Docker image.

## The CRAM-specific point that matters most

CRAM is *reference-based*: it cannot be decoded — and coverage cannot be
reproduced — without the exact reference it was compressed against. Each CRAM
`@SQ` header line carries an `M5` (MD5 of the sequence) tag.
`provenance.py::cram_reference_ids` extracts these so the supplied reference can
be verified against the one the CRAM was written with. That check is worth more
than hashing the 17 GB file (which the sidecar deliberately skips — it records
cheap size/mtime identity for big inputs and only sha256s the small `.crai`/`.fai`).

## Assumed CLI contract

`coverage.cwl` and the `Dockerfile` assume a console entrypoint that does not
exist yet. Adding it (`[project.scripts]` → `chromcov = "chromcov.cli:main"`)
would expose:

```
chromcov coverage \
  --cram data/COLO829T_TEST.cram \
  --reference data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa \
  --min-mapq 0 \
  --output coverage.tsv
```

...writing `coverage.tsv` and `coverage.tsv.provenance.json`. The CLI would call
the existing `chromcov.get_cov.calc_cov` per chromosome and
`provenance.build_provenance(...)` / `write_sidecar(...)` for the sidecar.

## Try the pieces now

```bash
# Provenance capture against this repo's own state (no CRAM needed):
python dev/reproducibility-sketch/provenance.py

# Build the pinned image:
docker build -t parklab-chromcov:0.1.0 -f dev/reproducibility-sketch/Dockerfile .

# Run via CWL and emit a CWLProv Research Object (RO-Crate) of the execution:
cwltool --provenance ro-crate/ \
  dev/reproducibility-sketch/coverage.cwl \
  dev/reproducibility-sketch/inputs.example.yml
```

## How this scales (for the write-up / Zoom)

Proportionate here; if this grew into a production pipeline: multi-step CWL
workflow on AWS (Batch/S3), GitHub Actions CI running the tool on a tiny subset
CRAM, and a heavier provenance/lineage layer (MLflow, or `redun`/`lamindb`) —
plus BioCompute Objects if it ever heads toward clinical/regulatory use.
