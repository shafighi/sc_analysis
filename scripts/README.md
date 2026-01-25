# Outlier Summary Pipeline

This document describes the pipeline for generating and visualizing outlier summaries from scAbsolute single-cell data.

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           OUTLIER SUMMARY PIPELINE                          │
└─────────────────────────────────────────────────────────────────────────────┘

Step 1: Prepare Sample CSV + Generate Outlier Summaries (per-sample RDS files)
         ↓
Step 2: Combine Summaries with Metadata (aggregated CSV)
         ↓
Step 3: Visualize Results (plots & figures)
```

---

## Step 1: Prepare Sample Manifest CSV & Generate Outlier Summaries

Create a CSV file listing the samples you want to process. Place it in `samples/` folder.

**Required columns:**
- `sample`: Sample identifier (e.g., `25393`, `23003`)

**Optional columns:**
- `category` or `Cell line`: Cell line name (e.g., `PEO1`, `HCT116`)
- `feature1`: First feature annotation (e.g., `HRD`, `HRP`)
- `feature2`: Second feature annotation (e.g., `Cisplatin sensitive`, `Cisplatin Resistant`)

**Example CSV (`samples/my_samples.csv`):**
```csv
sample,Cell line,feature1,feature2
23003,PEO1,HRD,Cisplatin sensitive
24077,PEO4,HRP,Cisplatin Resistant
24173,PEO14,HRD,Cisplatin Sensitive
```

See `samples/ALLSAMPLES_metadata.csv` for a complete example.

---

## Step 2: Generate Outlier Summaries

**Script:** `scripts/01_generate_outlier_summaries.R`

This script reads scAbsolute RDS objects and generates per-sample outlier summary files.

**Usage:**
```bash
Rscript scripts/01_generate_outlier_summaries.R <samples_csv> [obj_base] [out_base] [bin_size]
```

**Arguments:**
1. `samples_csv` - Path to sample manifest CSV
2. `obj_base` - (Optional) Directory containing input RDS files (default: `/Volumes/LenovoPS8/FI backup/sc_analysis/scAboslute-obj`)
3. `out_base` - (Optional) Output directory for results (default: `/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute`)
4. `bin_size` - (Optional) Bin size (default: `100`)

**Input:**
- Sample CSV manifest
- scAbsolute RDS objects at `obj_base/SLX-<sample>_<bin_size>.rds`

**Output (per sample):**
- `outlier_summary.rds` - Summary statistics
- `cellbased_outliers.rds` - Cell-level outlier data
- `<sample>_non_outlier.rds` - Filtered object without outliers
- `normals.rds` - Normal cells identified
- `heatmap_clustered.pdf` - Copy number heatmap

---

## Step 2: Combine Summaries with Metadata

**Script:** `scripts/02_combine_outlier_summaries.R`

This script aggregates individual `outlier_summary.rds` files into a single CSV with metadata.

**Usage:**
```bash
Rscript scripts/02_combine_outlier_summaries.R <samples_csv> <out_base> <obj_base> [bin_size]
```

**Arguments:**
1. `samples_csv` - Path to sample manifest CSV (with metadata columns)
2. `out_base` - Output directory for combined CSV
3. `obj_base` - Directory containing per-sample output folders (from Step 1)
4. `bin_size` - (Optional) Bin size used in processing (default: `100`)

**Output:**
- `outlier_summary_table_meta_combined.csv` - Combined summary with metadata
- `outlier_summary_table_meta.png` - Formatted table image (optional)

---

## Step 3: Visualize Results

**Script:** `scripts/03_visualize_summary.R`

A flexible, general-purpose visualizer that works with any summary CSV. Not specific to any sample type (FFPE, cell lines, etc.).

**Usage:**
```bash
Rscript scripts/03_visualize_summary.R <input_csv> <out_base> [group_col] [label_col]
```

**Arguments:**
1. `input_csv` - Path to combined outlier summary CSV (from Step 2)
2. `out_base` - Output directory for plots
3. `group_col` - (Optional) Column name for grouping/coloring (e.g., `feature1`, `category`)
4. `label_col` - (Optional) Column name for x-axis labels (e.g., `Sample`, `Cell line`)

**Output:**
- `<basename>_normalized.csv` - Normalized data with computed metrics
- `<basename>_quality_distribution.pdf` - Density plot of high-quality %
- `<basename>_quality_by_group.pdf` - Density plot colored by group
- `<basename>_quality_boxplot.pdf` - Box plot by group
- `<basename>_composition.pdf` - Stacked bar chart of cell composition
- `<basename>_summary_by_group.pdf` - Mean quality by group with error bars
- `<basename>_summary_stats.csv` - Summary statistics by group
- `<basename>_metrics_heatmap.pdf` - Heatmap overview of all metrics


## Quick Start Example

```bash
# Step 1: Generate outlier summaries for each sample
Rscript scripts/01_generate_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/scAboslute-obj" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute" \
    100

# Step 2: Combine all summaries into one CSV with metadata
Rscript scripts/02_combine_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute" \
    100

# Step 3: Generate visualization plots
Rscript scripts/03_visualize_summary.R \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output/outlier_summary_table_meta_combined.csv" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output" \
    feature1 \
    Sample
```

---

## Dependencies

Required R packages:
- `Biobase`
- `dplyr`
- `tidyr`
- `ggplot2`
- `knitr`
- `kableExtra` (optional, for PNG table rendering)
- `flextable` (optional, for PNG table rendering)
- `webshot2` (optional, for PNG export)
- `gridExtra` (for combined plots)

Install missing packages:
```r
install.packages(c("dplyr", "tidyr", "ggplot2", "knitr", "kableExtra", "flextable", "gridExtra"))
```

---

## Directory Structure

```
sc_analysis/
├── R/
│   ├── core.R                  # Core helper functions
│   ├── summary_helpers.R       # Summary computation helpers
│   └── visualization_helpers.R # Plotting helpers
├── samples/
│   ├── ALLSAMPLES_metadata.csv # Sample manifest with metadata
│   └── ...                     # Other sample lists
├── scripts/
│   ├── 01_generate_outlier_summaries.R
│   ├── 02_combine_outlier_summaries.R
│   ├── 03_visualize_outlier_summary.R
│   └── README.md               # This file
└── README.md                   # Project overview
```

---

## Troubleshooting

**"Object file not found" warning:**
- Check that `obj_base` path is correct and mounted
- Verify sample names match the RDS filenames

**"No outlier summaries found" error:**
- Ensure Step 2 completed successfully
- Check that `outlier_summary.rds` files exist in output folders

**PNG/PDF not generated:**
- Install optional packages: `kableExtra`, `webshot2`, `flextable`
- For `webshot2`, you may need to install Chrome/Chromium
