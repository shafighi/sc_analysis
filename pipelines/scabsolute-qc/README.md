# Post-scAbsolute QC pipeline

## What it does

The pipeline runs seven ordered stages:

1. Classify every cell as Replicating, low-RPC, metric outlier, Borderline, or
   PassedQC, and generate auditable per-cell copy-number profiles.
2. Combine per-sample QC summaries with manifest metadata.
3. Generate cross-sample QC composition plots.
4. Calculate per-cell and per-chromosome copy-number dropouts.
5. Plot the cross-sample dropout heatmap.
6. Plot per-sample dropout summaries and plate layouts.
7. Summarize total-read and RPC distributions.

## Run

From the repository root:

```bash
./run_scabsolute_qc.sh \
    <samples_csv> \
    <base_path> \
    [bin_size] \
    [qc_config] \
    [group_col] \
    [project_name]
```

Defaults:

- `bin_size`: `100`
- `qc_config`: `pipelines/scabsolute-qc/config/qc_params_default.csv`
- `group_col`: `Cell line`
- `project_name`: `all_samples`

Input objects must be named
`<base_path>/scAboslute-obj/SLX-<sample>_<bin_size>.rds`.

## QC order

Filters are sequential and categories are mutually exclusive:

1. Cells with unavailable HMM alpha are recorded separately.
2. Replicating/S-phase cells are identified and removed from metric QC.
3. Non-replicating cells below the RPC cutoff fail RPC QC.
4. MAPD and normalized Gini residuals are evaluated among RPC-passing cells.
5. HMM alpha is evaluated among RPC-passing cells.
6. Mild failures are assigned Borderline for manual review; stronger failures
   are outliers.
7. Remaining cells are PassedQC and may additionally be marked Normal.

The exact values used in a run are read from the selected config CSV and saved
with the outputs. Current maintained configs are:

- `config/qc_params_default.csv`
- `config/qc_params_relaxed.csv`

## Main outputs

Per sample, under `<base_path>/analysis_per_sample/SLX-<sample>_<bin_size>/`:

- `all_cells_qc.csv` - status, reason, and metrics for every cell.
- `qc_summary.csv` and `qc_params.csv` - counts and thresholds.
- `cells_passedqc.rds`, `cells_borderline.rds`, `cells_normal.rds`.
- `figures/cn_profiles_passedqc.pdf`.
- `figures/cn_profiles_outliers.pdf`.
- `figures/cn_profiles_borderline.pdf`.
- dropout and read-distribution tables/plots.

Across samples, under `<base_path>/results_<project_name>/`:

- combined QC CSV with run metadata.
- QC composition plots.
- dropout heatmap.

## Internal layout

- `run.sh` orchestrates the seven stages.
- `scripts/` contains only ordered pipeline stages.
- `R/` contains only helpers required by those stages.
- `config/` contains maintained QC configurations.
- `examples/` contains the manifest template.
