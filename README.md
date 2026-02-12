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
    100 \
    config/qc_params_default.csv

# Step 2: Combine all summaries into one CSV with metadata
Rscript scripts/02_combine_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute" \
    100

# Step 3: Generate visualization plots
Rscript scripts/03_visualize_summary.R \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output/{output}_{config}_{date}.csv" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output" \
    feature1 \
    Sample

# Step 4: Generate per-sample dropout summaries (requires step 1 output)
Rscript scripts/04_generate_dropout_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/scAboslute-obj" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute" \
    100

# Step 5: Visualize dropout heatmap (requires step 4 output)
Rscript scripts/05_visualize_dropout_heatmap.R \
    samples/ALLSAMPLES_metadata.csv \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/post-scAbsolute" \
    "/Volumes/LenovoPS8/FI backup/sc_analysis/all_samples_23july2024/output" \
    100 \
    "Cell line"
```

See [scripts/README.md](scripts/README.md) for detailed pipeline documentation.

## Output

**Output file** (named by output folder + config + date):
- `{output_folder}_{config}_{YYYYMMDD}.csv` - Combined QC summary with metadata header

Same config + same date = overwrites previous file (no redundant copies).

**Output file format:**
```
# ==============================================================================
# QC Summary Report
# ==============================================================================
# Run Date: 20240125
# Samples CSV: ALLSAMPLES_metadata.csv
# Samples Processed: 22
# Bin Size: 100
#
# QC Parameters:
#   QC_RPC_cutoff: 25
#   QC_MAPD_cutoff: 2
#   QC_Gini_cutoff: 2
#   QC_Alpha_cutoff: 1.5
#   ...
# ==============================================================================
"Sample","Cell line","feature1",...
"23003","PEO1","HRD",...
```

## QC Summary Columns

| Column | Description |
|--------|-------------|
| `Sample` | Sample identifier |
| `Cell line` | Cell line name |
| `feature1` | First annotation (e.g., HRD/HRP) |
| `feature2` | Second annotation (e.g., treatment sensitivity) |
| `post-scAbsolute` | Cells that passed scAbsolute processing (not all sequenced cells) |
| `Replicating` | Cells in S-phase (identified by read depth variation) |
| `Outliers(RPC)` | Non-replicating cells with RPC < 25 (filtered first) |
| `Outliers(Alpha/Mapd/Gini,post-RPC)` | Cells that passed RPC but failed Alpha/Mapd/Gini QC |
| `PassedQC(incl.Normal)` | Cells passing all QC filters (includes Normal) |
| `Normal` | Diploid cells where >95% of bins have CN=2 |
| `PassedQC-Normal` | Aberrant cells with CNVs that passed QC |
| `(PassedQC+Replicating)/post-scAbsolute%` | Percentage of usable cells |
| `PassedQC/post-scAbsolute%` | Percentage of cells that passed all QC |
| `(PassedQC-Normal)/post-scAbsolute%` | Percentage of aberrant (CNV) cells |

**Note:** QC parameters (RPC cutoff, MAPD cutoff, etc.) are stored in the file header as comments, not as columns.

## Cell Flow

```
post-scAbsolute (cells passed scAbsolute)
├── Replicating (S-phase cells, set aside)
└── Non-replicating
    ├── Outliers(RPC) ──────────────────────► removed (RPC < 25)
    └── Passed RPC
        ├── Outliers(Alpha/Mapd/Gini) ──────► removed (QC failed)
        └── PassedQC(incl.Normal)
            ├── Normal (diploid, >95% CN=2)
            └── PassedQC-Normal (aberrant, has CNVs)
```

**Formula:** `post-scAbsolute` = `Replicating` + `Outliers(RPC)` + `Outliers(Alpha/Mapd/Gini)` + `PassedQC(incl.Normal)`

## QC Pipeline Explanation

The outlier detection pipeline applies **sequential filtering** on non-replicating cells:

1. **Replicating cells are identified first** based on read depth variation patterns indicating S-phase. These cells are set aside (not removed, but tracked separately).

2. **RPC (Reads Per Cell) filtering is applied first** to non-replicating cells. Cells with RPC < 25 are marked as `Outliers(RPC)` and removed from further QC steps. This ensures sufficient sequencing depth for reliable copy number calls.

3. **Alpha/MAPD/Gini filtering is applied second**, only to cells that passed RPC:
   - **MAPD** (Median Absolute Pairwise Difference): measures noise between adjacent bins
   - **Gini coefficient**: measures inequality in read distribution
   - **HMM Alpha**: measures segmentation confidence

   Cells failing any of these metrics are marked as `Outliers(Alpha/Mapd/Gini,post-RPC)`.

4. **PassedQC** cells are those that passed ALL filters (RPC + Alpha/MAPD/Gini). This includes both aberrant cells (with CNVs) and normal diploid cells.

5. **Normal cells** are a subset of PassedQC where >95% of genomic bins have copy number = 2 (diploid).

6. **PassedQC-Normal** = aberrant cells that passed QC but have copy number variations.

**Key point:** The outlier categories are **non-overlapping** because filtering is sequential. A cell marked as `Outliers(RPC)` never gets evaluated for Alpha/MAPD/Gini, so there's no double-counting.

## QC Configuration

QC parameters are stored in CSV config files in the `config/` directory:

| Config File | Description |
|-------------|-------------|
| `config/qc_params_default.csv` | Standard thresholds (recommended) |
| `config/qc_params_relaxed.csv` | More permissive thresholds for exploratory analysis |

### Default QC Thresholds

| Parameter | Default | Relaxed | Description |
|-----------|---------|---------|-------------|
| `rpc_cutoff` | 25 | 15 | Reads Per Cell minimum |
| `mapd_cutoff` | 2.0 | 1.5 | MAPD maximum |
| `gini_norm_cutoff` | 2.0 | 1.5 | Gini coefficient maximum |
| `alpha_cutoff` | 1.5 | 2.0 | HMM Alpha (SD multiplier) |
| `alpha_hard_cutoff` | 0.05 | 0.05 | HMM Alpha absolute threshold |
| `normal_threshold` | 95% | 90% | % bins with CN=2 for normal cells |

### Using Different QC Parameters

```bash
# Use default (strict) parameters
Rscript scripts/01_generate_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/scAbsolute-obj \
    /path/to/output \
    100 \
    config/qc_params_default.csv

# Use relaxed parameters for exploratory analysis
Rscript scripts/01_generate_outlier_summaries.R \
    samples/ALLSAMPLES_metadata.csv \
    /path/to/scAbsolute-obj \
    /path/to/output \
    100 \
    config/qc_params_relaxed.csv
```

### Creating Custom QC Parameters

Copy and modify an existing config file:

```bash
cp config/qc_params_default.csv config/qc_params_custom.csv
# Edit config/qc_params_custom.csv with your values
```

Config file format (CSV):
```csv
parameter,value,description
rpc_cutoff,25,Reads Per Cell minimum
mapd_cutoff,2,MAPD maximum
gini_norm_cutoff,2,Gini coefficient maximum
alpha_cutoff,1.5,HMM Alpha cutoff (SD multiplier)
...
```

The QC parameters used are automatically saved in `qc_criteria_used.csv` per sample and in the header of the combined output CSV.

## Dropout Analysis (Steps 4-5)

Steps 4 and 5 analyze **copy number dropouts** — genomic bins where the segmented copy number is 0, indicating potential technical artifacts or biological deletions.

### Step 4: Generate Dropout Summaries

For each sample, reads the scAbsolute object and `cellbased_outliers.rds` from step 1, then:
- Counts dropout bins (segVal == 0) per cell per chromosome
- Assigns QC status to each cell (pass, replicating, rpc_outlier, etc.)
- Saves per-sample files: `cn_binned.rds`, `cell_chr_dropout.rds`, `cell_dropout_status.csv`
- Writes a combined `dropout_summary_combined.csv` with per-chromosome dropout frequencies normalized by cell count

### Step 5: Visualize Dropout Heatmap

Generates a publication-ready heatmap (PDF) showing dropout frequency per cell across chromosomes and samples:
- Compact layout with numeric annotations in each tile
- Viridis color scale, faceted by sample category
- Adaptive text color (white on dark, black on light tiles)
- Auto-sized dimensions based on number of samples

Both steps can be run independently or as part of the full pipeline.

## Additional Analysis Scripts

### Cluster Classification Analysis

Machine learning classification and clustering for comparing cell populations across conditions (e.g., HRP vs IHR samples).

**Features:**
- LASSO regression for identifying discriminative genomic positions
- t-SNE dimensionality reduction for visualization
- PCA + hierarchical clustering for sample grouping

```bash
Rscript scripts/cluster_classification_analysis.R \
    <input_dir> \
    <output_dir> \
    <bin_size> \
    <sample_ids> \
    <labels>

# Example: Compare IHR (label=1) vs HRP (label=0) samples
Rscript scripts/cluster_classification_analysis.R \
    /path/to/scAbsolute-obj \
    /path/to/output \
    100 \
    "24490,24489,24174,24175" \
    "1,0,1,0"
```

**Output files:**
- `lasso_important_positions.csv` - Genomic positions distinguishing groups
- `tsne_plot.pdf/png` - t-SNE visualization
- `dendrogram.pdf` - Hierarchical clustering tree
- `pca_clusters.pdf/png` - PCA with cluster assignments
- `cluster_results.csv` - Per-cell cluster assignments
