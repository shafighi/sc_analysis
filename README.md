# sc_analysis

Single-cell copy number analysis pipeline for processing scAbsolute and scUnique results.

## Quick Start

Run the complete pipeline with a single command:

```bash
./run_pipeline.sh
```

Or with custom parameters:

```bash
./run_pipeline.sh [samples_csv] [base_path] [bin_size]

# Examples:
./run_pipeline.sh                                                    # use defaults
./run_pipeline.sh samples/my_samples.csv                             # custom sample list
./run_pipeline.sh samples/my.csv "/Volumes/MyDrive/data" 50          # custom path and bin size
```

**Default values:**
- `samples_csv`: `samples/ALLSAMPLES_metadata.csv`
- `base_path`: `/Volumes/LenovoPS8/FI backup/sc_analysis`
- `bin_size`: `100`

## Individual Steps

You can also run each step separately:

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

See [scripts/README.md](scripts/README.md) for detailed pipeline documentation.

## Output Column Definitions

### Raw Columns (from `outlier_summary_table_meta_combined.csv`)

| Column | Description |
|--------|-------------|
| `Processed Cells` | Total number of cells sequenced in the sample |
| `Replicating` | Cells detected as replicating (S-phase) based on read depth variation |
| `Replicating & RPC` | Replicating cells that also fall below the RPC threshold |
| `RPC Outliers` | Non-replicating cells with reads per cell < 25 (low coverage) |
| `Alpha/Mapd/Gini` | Cells failing any QC metric (HMM alpha, MAPD, or Gini coefficient) |
| `Good Quality Cells` | Cells passing all QC filters (non-replicating, sufficient RPC, good alpha/MAPD/Gini) |
| `Normal Cells` | Diploid cells where >95% of bins have copy number = 2 |
| `Sample` | Sample identifier |
| `Cell line` | Cell line name from metadata |
| `feature1` | First annotation (e.g., HRD/HRP status) |
| `feature2` | Second annotation (e.g., treatment sensitivity) |

### Derived Columns (in `*_normalized.csv`)

| Column | Formula | Description |
|--------|---------|-------------|
| `Sequenced` | = `Processed Cells` | Renamed for clarity |
| `Good_Quality` | = `Good Quality Cells` | Renamed for clarity |
| `Normal` | = `Normal Cells` | Renamed for clarity |
| `Outliers` | = `Alpha/Mapd/Gini` | Renamed for clarity |
| `Filtered` | = `Good_Quality` - `Normal` | Good quality cells excluding normal diploid cells |
| `High_Quality_Pct` | = (`Good_Quality` + `Replicating`) × 100 / `Sequenced` | Percentage of cells that are usable (good quality + replicating) |
| `Pass_Rate` | = (`Sequenced` - `Outliers`) × 100 / `Sequenced` | Percentage of cells not flagged as QC outliers |
| `Label` | from `label_col` parameter | X-axis label for plots (default: Sample) |
| `Group` | from `group_col` parameter | Grouping variable for comparisons (default: feature1) |
