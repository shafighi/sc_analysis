# Repository guidance

## Scope

This repository has exactly two maintained pipelines:

1. `pipelines/scabsolute-qc/` - post-scAbsolute cell QC and reporting.
2. `pipelines/post-scunique/` - post-scUnique tree/CN and private-event plots.

Everything under `archive/` is historical reference material. Do not source,
extend, or treat archived code as a dependency of the maintained pipelines.

## Entry points

- `run_scabsolute_qc.sh`
- `run_post_scunique.sh`

Keep these root launchers thin. Pipeline logic belongs in the corresponding
`pipelines/<name>/` directory.

## scAbsolute QC pipeline

- `run.sh` orchestrates seven ordered scripts.
- `scripts/` contains only numbered maintained stages.
- `R/` contains only helpers used by those stages.
- `config/` contains maintained threshold sets.
- Scripts execute with `pipelines/scabsolute-qc/` as their working directory,
  so internal `source("R/...")` paths are pipeline-relative.

QC order is replication, RPC, MAPD/Gini/alpha, Borderline review, then
PassedQC/Normal classification. Preserve mutually exclusive final categories
and explicit per-cell failure reasons.

## post-scUnique pipeline

The input contract is a sample result directory containing matching
`finalCN.RDS`, `tree.RDS`, and `df_pass_post.RDS` files. The pipeline must keep
full cell names in CSV output even when plot labels are shortened.

`freq == 1` means a validated locus is carried by one cell in the analyzed
cohort. A cell can carry several different private events.

## Validation

For structural changes:

```bash
bash -n run_scabsolute_qc.sh run_post_scunique.sh \
  pipelines/scabsolute-qc/run.sh pipelines/post-scunique/run.sh

Rscript --vanilla -e 'parse(file="path/to/script.R")'
```

When changing post-scUnique plotting, run it on a completed sample and render
the PDFs to PNG for visual inspection. When changing scAbsolute orchestration,
verify all seven stage calls and run at least one real sample through stage 1.

## Archive policy

Move superseded or one-off material into the appropriate `archive/` category
instead of mixing it into either active pipeline. Preserve archived files unless
the user explicitly asks to delete them.
