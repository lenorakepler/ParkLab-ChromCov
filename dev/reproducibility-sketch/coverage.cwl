#!/usr/bin/env cwl-runner
# CWL CommandLineTool wrapping the chromcov coverage step.
#
# Run with the CWL reference runner and get a CWLProv Research Object (RO-Crate)
# of the whole execution for free:
#
#     cwltool --provenance ro-crate/ coverage.cwl inputs.example.yml
#
# NOTE: this wraps an assumed `chromcov coverage` CLI (see README) that does not
# exist yet — it defines the contract the CLI should expose.
cwlVersion: v1.2
class: CommandLineTool
label: Per-chromosome average coverage from a CRAM

# The Docker image pins samtools/htslib + pysam + the chromcov code itself, so
# the environment is reproducible independent of the host.
requirements:
  DockerRequirement:
    dockerPull: ghcr.io/lenorakepler/parklab-chromcov:0.1.0

baseCommand: [chromcov, coverage]

inputs:
  cram:
    type: File
    # CRAM decoding needs the index; CWL stages secondaryFiles alongside it.
    secondaryFiles:
      - pattern: .crai
        required: true
    inputBinding:
      prefix: --cram
  reference:
    type: File
    secondaryFiles:
      - pattern: .fai
        required: true
    inputBinding:
      prefix: --reference
  min_mapq:
    type: int
    default: 0
    inputBinding:
      prefix: --min-mapq
  output_name:
    type: string
    default: coverage.tsv
    inputBinding:
      prefix: --output

outputs:
  coverage_table:
    type: File
    outputBinding:
      glob: $(inputs.output_name)
  provenance:
    # Written by provenance.py next to the table: <output>.provenance.json
    type: File
    outputBinding:
      glob: $(inputs.output_name).provenance.json
