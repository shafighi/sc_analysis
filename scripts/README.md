# Outlier Summary Pipeline

This document describes the pipeline for generating and visualizing outlier summaries from scAbsolute single-cell data.

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           OUTLIER SUMMARY PIPELINE                          │
└─────────────────────────────────────────────────────────────────────────────┘

Step 1: Prepare Sample CSV
         ↓
Step 2: Generate Outlier Summaries (per-sample RDS files)
         ↓
Step 3: Combine Summaries with Metadata (aggregated CSV)
         ↓
Step 4: Visualize Results (plots & figures)
```

---

## Step 1: Prepare Sample Manifest CSV

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
Rscript scripts/01_generate_outlier_summaries.R samples/my_samples.csv
```

**Input:**
- Sample CSV manifest
- scAbsolute RDS objects at `obj_base/SLX-<sample>_<bin_size>.rds`

**Output (per sample):**
- `outlier_summary.rds` - Summary statistics
- `cellbased_outliers.rds` - Cell-level outlier data
- `<sample>_non_outlier.rds` - Filtered object without outliers
- `normals.rds` - Normal cells identified
- `heatmap_clustered.pdf` - Copy number heatmap

**Configuration:**
Edit the script to adjust paths:
```r
obj_base <- "/Volumes/Fl/sc_analysis/scAboslute-obj"  # Input RDS location
out_base <- "/Volumes/Fl/sc_analysis/post-scAbsolute" # Output directory
bin_size <- "100"                                      # Bin size
```

---

## Step 3: Combine Summaries with Metadata

**Script:** `scripts/02_combine_outlier_summaries.R`

This script aggregates individual `outlier_summary.rds` files into a single CSV with metadata.

**Usage:**
```bash
Rscript scripts/02_combine_outlier_summaries.R samples/my_samples.csv /path/to/output /path/to/post-scAbsolute 100
```

**Arguments:**
1. `samples_csv` - Path to sample manifest CSV (with metadata columns)
2. `out_base` - Output directory for combined CSV
3. `obj_base` - Directory containing per-sample output folders
4. `bin_size` - Bin size used in processing (default: 100)

**Output:**
- `outlier_summary_table_meta_combined.csv` - Combined summary with metadata
- `outlier_summary_table_meta.png` - Formatted table image (optional)

---

## Step 4: Visualize Results

**Script:** `scripts/03_visualize_outlier_summary.R`

This script generates visualization plots from the combined CSV.

**Usage:**
```bash
Rscript scripts/03_visualize_outlier_summary.R /path/to/outlier_summary_table_meta_combined.csv /path/to/ffpe.csv /path/to/output
```

**Arguments:**
1. `combined_csv` - Path to combined outlier summary CSV (from Step 3)
2. `ffpe_csv` - Path to FFPE samples CSV (optional, pass empty string to skip)
3. `out_base` - Output directory for plots

**Output:**
- `Good_Quality_Distribution.pdf` - Density plot of high-quality cell percentage
- `cell_composition_plot.pdf` - Stacked bar chart of cell composition
- `Cellines_95_2.csv` - Normalized data with computed percentages
- For FFPE samples (if provided):
  - `Good_Quality_Distribution_FFPE.pdf`
  - `cell_composition_plot_FFPE.pdf`
  - `Good_Quality_Distribution_Combined.pdf`

---

## Complete Example Workflow

```bash
# 1. Prepare your sample manifest
cp samples/ALLSAMPLES_metadata.csv samples/my_analysis.csv
# Edit my_analysis.csv to include your samples

# 2. Generate outlier summaries for each sample
Rscript scripts/01_generate_outlier_summaries.R samples/my_analysis.csv

# 3. Combine all summaries into one table with metadata
Rscript scripts/02_combine_outlier_summaries.R \
    samples/my_analysis.csv \
    /Volumes/Fl/sc_analysis/results \
    /Volumes/Fl/sc_analysis/post-scAbsolute \
    100

# 4. Generate visualization plots
Rscript scripts/03_visualize_outlier_summary.R \
    /Volumes/Fl/sc_analysis/results/outlier_summary_table_meta_combined.csv \
    "" \
    /Volumes/Fl/sc_analysis/results
```

---

## File Naming Convention

| Old Name | New Name | Purpose |
|----------|----------|---------|
| `generate_summary_of_outliers_scAbsolute.R` | `01_generate_outlier_summaries.R` | Generate per-sample outlier RDS files |
| `visualize_summary_outlier_metadata_table.R` | `02_combine_outlier_summaries.R` | Combine RDS files into CSV with metadata |
| `visualize_summary.R` | `03_visualize_outlier_summary.R` | Create plots from combined CSV |
| `visualize_summary_outlier_table.R` | (deprecated) | Simpler version without full metadata |

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
