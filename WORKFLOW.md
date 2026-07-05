# chromcov workflow

Coverage is computed once and reused. The pipeline is **three derivation levels**,
each writing real, content-addressed output files, so a later call checks what
already exists and picks up where it left off (the "trickle-down" reuse).

## Data flow + reuse

```mermaid
flowchart TD
    CRAM[("CRAM + .crai + reference")]:::in
    BEDS[("SMaHT strata BEDs")]:::in
    FETCH["chromcov fetch strata"] -.->|downloads| BEDS

    CRAM --> PRE{"preflight<br/>sorted · indexed · reference-M5 match"}

    PRE --> DECIDE{"per-base track exists?<br/>key = inputs + read-filter params"}
    DECIDE -->|"miss"| CALC["calc_cov(per_base=True)<br/>→ per-base depth vector"]
    DECIDE -->|"hit — reuse"| LOAD["load track<br/>→ reconstruct depth vector"]

    CALC -->|"--per-base writes it"| TRACK[("Level 1 · per-base depth tracks<br/>out/perbase/&lt;coverage-key&gt;/chrN.per-base.bedgraph.gz<br/>+ coverage.json sidecar")]:::out
    TRACK -.->|reused next run| DECIDE
    TRACK -.->|"interoperable output"| EXT["bedtools · IGV · bigWig"]

    CALC --> RED["ChromDepth reductions<br/>histogram · windows · strata masks"]
    LOAD --> RED
    BEDS --> RED

    RED --> ANALYSIS[("Level 2 · analysis run<br/>out/analysis/&lt;analysis-key&gt;/<br/>stats · windows · strata · plots + run.json<br/>key = coverage-key + analysis params")]:::out

    ANALYSIS -.->|deferred| COLLATE["Level 3 · collate/compare runs<br/>(stratified vs not, ...)"]:::todo

    classDef in fill:#e8f0fe,stroke:#4c72b0,color:#000;
    classDef out fill:#e6f4ea,stroke:#34a853,color:#000;
    classDef todo fill:#fef7e0,stroke:#f9ab00,stroke-dasharray:4 3,color:#000;
```

**Why the keys differ per level.** The expensive step (per-base depth) depends only
on the *inputs + read-filter params* — so its key ignores window size, strata, or
copy-number settings. Changing `--window` or adding `--strata` therefore reuses the
same Level-1 tracks and only re-runs the cheap Level-2 reductions. That's what makes
"stratified vs unstratified" a fast comparison: two Level-2 runs off one Level-1 key.

**`coverage` (the deliverable) is the shortcut path:** it takes `per_base=False`,
sums aligned bases straight to `mean = bases / length`, and never builds the vector
or the tracks — fast, for the headline table (and the CWL contract).

## Orchestration layers

The same steps are usable three ways: call the library, run the CLI, or let an
engine drive the CLI. The per-base tracks are the file substrate an engine checks.

```mermaid
flowchart LR
    LIB["chromcov package<br/>CoverageConfig · backends · ChromDepth<br/>Strata · CoverageAnalysis · PerBaseStore"] --> CLI["chromcov CLI<br/>coverage · perbase · analyze · fetch · collate"]
    CLI --> SNAKE["Snakefile<br/>rule per chromosome<br/>scatter tracks → gather analysis<br/>(parallelism · resume · --dry-run)"]
    CLI --> CWLW["coverage.cwl + Docker<br/>portable tool wrapper"]

    classDef n fill:#f5f5f5,stroke:#888,color:#000;
    class LIB,CLI,SNAKE,CWLW n;
```

- **Built-in reuse** (the decision diamond above) makes `chromcov` self-contained —
  no engine needed to get incremental behavior.
- **Snakemake** scatters one `chromcov perbase --chrom N` per chromosome (free
  parallelism + resume), then gathers into one `chromcov analyze` that reuses the
  tracks. It imports `chromcov` to compute the deterministic keys, so its rule
  outputs land on the exact same content-addressed paths.
- **CWL** wraps `chromcov coverage` as a portable, Dockerized tool for the SMaHT
  pipeline world (the JD's stated target).
