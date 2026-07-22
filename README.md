# sc_analysis

Two focused pipelines for single-cell copy-number analysis:

| Pipeline | Purpose | Entry point |
|---|---|---|
| post-scAbsolute QC | Classify cells, create auditable QC reports, and summarize dropout/read distributions | `./run_scabsolute_qc.sh` |
| post-scUnique visualization | Plot the evolutionary tree with final copy number and summarize private (`freq == 1`) events per cell | `./run_post_scunique.sh` |

Historical analyses, old plotting scripts, and previous sample manifests are
kept under [`archive/`](archive/README.md). They are preserved for reference but
are not part of either maintained pipeline.

## Repository layout

```text
sc_analysis/
├── run_scabsolute_qc.sh
├── run_post_scunique.sh
├── environment.yml
├── pipelines/
│   ├── scabsolute-qc/
│   │   ├── run.sh
│   │   ├── R/
│   │   ├── scripts/
│   │   ├── config/
│   │   └── examples/
│   └── post-scunique/
│       ├── run.sh
│       └── scripts/
└── archive/
```

## Installation

Create the shared Conda environment:

```bash
conda env create -f environment.yml
conda activate sc_analysis
```

Set `RSCRIPT_BIN` if `Rscript` is not on `PATH`:

```bash
export RSCRIPT_BIN=/path/to/Rscript
```

## Pipeline 1: post-scAbsolute QC

This pipeline consumes scAbsolute `QDNAseqCopyNumbers` RDS objects. It performs
cell-cycle classification, RPC filtering, MAPD/Gini/alpha QC, Borderline
classification, manual-review profile reports, dropout summaries, and
cross-sample reporting.

```bash
./run_scabsolute_qc.sh \
    path/to/samples.csv \
    path/to/project_data \
    500
```

The project data directory must contain:

```text
project_data/
└── scAboslute-obj/
    └── SLX-<sample>_<bin_size>.rds
```

Results are written to:

- `project_data/analysis_per_sample/` - cell classifications, QC tables, and
  per-sample PDFs.
- `project_data/results_<project_name>/` - combined tables and cross-sample
  figures.

See [`pipelines/scabsolute-qc/README.md`](pipelines/scabsolute-qc/README.md) for
arguments, QC order, thresholds, and outputs.

## Pipeline 2: post-scUnique visualization

This pipeline consumes one completed scUnique sample directory and produces a
tree-aligned final-copy-number heatmap plus the distribution of private events.
Here, `freq == 1` means an event locus was observed in exactly one cell in the
analyzed scUnique cohort; a cell may carry several such private events.

```bash
./run_post_scunique.sh \
    path/to/scunique/results/SLX-27548_500 \
    path/to/output
```

Required files:

- `<prefix>.finalCN.RDS`
- `<prefix>.tree.RDS`
- `<prefix>.df_pass_post.RDS`

See [`pipelines/post-scunique/README.md`](pipelines/post-scunique/README.md) for
all generated figures and tables.

## Sample manifest

The scAbsolute pipeline requires a CSV with a `sample` column. Optional metadata
columns are used for grouped summaries:

```csv
sample,Cell line,feature1,feature2
27548,Example line,Group A,Condition A
```

Start from
[`pipelines/scabsolute-qc/examples/sample_manifest.csv`](pipelines/scabsolute-qc/examples/sample_manifest.csv).
