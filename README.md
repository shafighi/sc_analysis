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

### Normalized Columns (in `*_normalized.csv`)

| Column | Description |
|--------|-------------|
| `Sample` | Sample identifier |
| `Cell line` | Cell line name |
| `feature1` | First annotation (e.g., HRD/HRP) |
| `feature2` | Second annotation (e.g., treatment sensitivity) |
| `post-scAbsolute` | Total cells sequenced (= `Processed Cells`) |
| `Replicating` | Cells in S-phase |
| `Outliers(RPC)` | Non-replicating cells with low read count (RPC < 25). Filtered FIRST, no overlap with Alpha/Mapd/Gini outliers |
| `Outliers(Alpha/Mapd/Gini,post-RPC)` | Cells that passed RPC but failed Alpha/Mapd/Gini QC. No overlap with RPC outliers |
| `PassedQC(incl.Normal)` | Cells passing all QC filters. Includes Normal cells as a subset |
| `Normal` | Diploid cells where >95% of bins have copy number = 2 (subset of PassedQC) |
| `PassedQC-Normal` | Aberrant cells with CNVs that passed QC |
| `(PassedQC+Replicating)/post-scAbsolute%` | Percentage of usable cells |
| `PassedQC/post-scAbsolute%` | Percentage of cells that passed all QC |
| `(PassedQC-Normal)/post-scAbsolute%` | Percentage of aberrant (CNV) cells |

**Cell Flow:**
```
post-scAbsolute
├── Replicating (S-phase cells, set aside)
└── Non-replicating
    ├── Outliers(RPC) ──────────────────────► removed (low read count)
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

   Cells failing any of these metrics are marked as `Outliers(Alpha/Mapd/Gini)`.

4. **PassedQC** cells are those that passed ALL filters (RPC + Alpha/MAPD/Gini). This includes both aberrant cells (with CNVs) and normal diploid cells.

5. **Normal cells** are a subset of PassedQC where >95% of genomic bins have copy number = 2 (diploid).

6. **PassedQC-Normal** = aberrant cells that passed QC but have copy number variations.

**Key point:** The outlier categories are **non-overlapping** because filtering is sequential. A cell marked as `Outliers(RPC)` never gets evaluated for Alpha/MAPD/Gini, so there's no double-counting.
