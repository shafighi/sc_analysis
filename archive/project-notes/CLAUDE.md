# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Single-cell copy number analysis pipeline for processing results from scAbsolute and scUnique bioinformatics tools. Focuses on identifying outlier cells in single-cell DNA sequencing data and generating visualizations of copy number variations.

## Architecture

### Core Components

- **R/** - Shared helper modules
  - `core.R` - Main QC functions: `qc_gini()`, `qc_alpha()`, `qc_mapd()`, `predict_replicating()`, `get_summary_of_outliers()`
  - `summary_helpers.R` - Outlier summarization
  - `visualization_helpers.R` - `plotCopynumberHeatmap()` and plotting utilities
  - `dropout_helpers.R` - Dropout analysis: `get_dropout_by_chromosome()`, `assign_cell_status()`

- **scripts/** - Analysis pipeline (numbered scripts are the main workflow)
  - `01_generate_outlier_summaries.R` - Process individual scAbsolute RDS objects, apply QC filters
  - `02_combine_outlier_summaries.R` - Aggregate summaries with metadata into combined CSV
  - `03_visualize_summary.R` - Generate publication-ready plots
  - `04_generate_dropout_summaries.R` - Analyze copy number dropouts per cell and chromosome
  - `05_visualize_dropout_heatmap.R` - Cross-sample dropout frequency heatmap
  - `06_visualize_sample_dropouts.R` - Per-sample dropout barplots and cell x chromosome heatmaps

- **samples/** - Sample manifest CSVs with metadata (sample ID, cell line, features)

### Directory Structure

```
{base_path}/
├── scAboslute-obj/              # Input (scAbsolute RDS objects)
├── analysis_per_sample/         # Per-sample outputs (steps 1, 4, 6)
│   ├── SLX-{sample}_{bin}/
│   │   ├── cellbased_outliers.rds
│   │   ├── cn_binned.rds
│   │   ├── figures/
│   │   └── ...
│   └── dropout_summary_combined.csv
└── results_{project_name}/      # Cross-sample outputs (steps 2, 3, 5)
    ├── output_{config}_{date}.csv
    ├── *_composition.pdf
    └── dropout_heatmap_{date}.pdf
```

### Data Flow

Input: scAbsolute RDS objects (`SLX-{sample}_{bin_size}.rds`) containing QDNAseq copy number objects
Output: Filtered objects, quality summaries, heatmap PDFs, aggregated CSV, visualization plots

## Running the Pipeline

```bash
# Full pipeline via run_pipeline.sh (recommended)
# Usage: ./run_pipeline.sh [samples_csv] [base_path] [bin_size] [qc_config] [group_col] [project_name]
./run_pipeline.sh samples/ALLSAMPLES_metadata.csv /path/to/data 100
./run_pipeline.sh samples/PEO_samples.csv /path/to/data 100 config/qc_params_default.csv "Cell line" peo

# Individual steps (for standalone use):

# Step 1: Generate per-sample outlier summaries
Rscript scripts/01_generate_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/scAbsolute-obj \
    /path/to/analysis_per_sample \
    100

# Step 2: Combine summaries with metadata
Rscript scripts/02_combine_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/results_all_samples \
    /path/to/analysis_per_sample \
    100

# Step 3: Visualize combined summary
Rscript scripts/03_visualize_summary.R \
    /path/to/results_all_samples/output_*.csv \
    /path/to/results_all_samples \
    feature1 \
    Sample

# Step 4: Generate per-sample dropout summaries (requires step 1 output)
Rscript scripts/04_generate_dropout_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/scAbsolute-obj \
    /path/to/analysis_per_sample \
    100

# Step 5: Visualize dropout heatmap (requires step 4 output)
Rscript scripts/05_visualize_dropout_heatmap.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/analysis_per_sample \
    /path/to/results_all_samples \
    100 \
    "Cell line"

# Step 6: Per-sample dropout plots (requires step 4 output)
Rscript scripts/06_visualize_sample_dropouts.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/analysis_per_sample \
    /path/to/analysis_per_sample \
    100
```

## QC Filtering Order

Filters are applied sequentially. Each step operates on cells that passed the previous step:

1. **NA Alpha removal** — Cells with missing `hmm.alpha` values are excluded from QC (tracked separately)
2. **Replicating cell detection** — S-phase/cycling cells identified via `predict_replicating()` and removed from further QC
3. **RPC filter** — Non-replicating cells with `rpc < rpc_cutoff` are flagged as `Outliers(RPC)` and excluded
4. **MAPD outlier detection** — Robust regression of MAPD ~ 1/RPC; residuals > cutoff flagged (`qc_mapd()`)
5. **Gini outlier detection** — Robust regression of Gini ~ RPC; residuals > cutoff flagged (`qc_gini_norm()`)
6. **Alpha outlier detection** — Log-normal distribution; cells beyond cutoff SD or hard threshold flagged (`qc_alpha()`)
7. **Combined outlier** — Cells failing any of MAPD/Gini/Alpha (logical OR) are `Outliers(Alpha/Mapd/Gini,post-RPC)`
8. **Normal cell detection** — Among PassedQC cells, those with >95% bins at CN=2 are flagged as `Normal`

Summary categories (mutually exclusive, sum to `post-scAbsolute`):
- `Replicating` — cycling cells (includes `Replicating(low-RPC)` subset)
- `Outliers(RPC)` — non-replicating cells with insufficient reads
- `Outliers(Alpha/Mapd/Gini,post-RPC)` — quality outliers detected after RPC filtering
- `PassedQC(incl.Normal)` — cells passing all filters (includes `Normal` subset)

## QC Metrics and Default Thresholds

- **RPC** (Reads Per Cell): >= 25
- **MAPD** (Median Absolute Pairwise Difference): <= 2.0
- **Gini Coefficient**: normalized <= 2.0
- **HMM Alpha**: <= 1.5 (or hard cutoff 0.05)

## Key Dependencies

Biobase, QDNAseq, ComplexHeatmap, circlize, viridis, dplyr, tidyr, ggplot2, ggbeeswarm, ggpubr, robustbase, gridExtra

## Sample Manifest Format

CSV with required `sample` column and optional metadata columns (`Cell line`, `feature1`, `feature2`):
```csv
sample,Cell line,feature1,feature2
23003,PEO1,HRD,Cisplatin sensitive
```
