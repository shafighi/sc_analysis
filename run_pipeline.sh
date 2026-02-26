#!/bin/bash
# ==============================================================================
# Run the complete analysis pipeline (Steps 1-7)
# ==============================================================================
#
# Usage: ./run_pipeline.sh [samples_csv] [base_path] [bin_size] [qc_config] [group_col] [project_name]
#
# Steps:
#   1 - Generate per-sample outlier summaries
#   2 - Combine summaries with metadata
#   3 - Visualize combined summary
#   4 - Generate per-sample dropout summaries
#   5 - Cross-sample dropout heatmap
#   6 - Per-sample dropout barplots and heatmaps
#   7 - Generate read distribution summaries
#
# Arguments:
#   samples_csv - CSV file with sample IDs (default: samples/ALLSAMPLES_metadata.csv)
#   base_path   - Base path to data directories (default: /Volumes/LenovoPS8/FI backup/sc_analysis)
#   bin_size    - Bin size for analysis (default: 100)
#   qc_config   - QC parameters config file (default: config/qc_params_default.csv)
#   group_col    - Column to group/facet dropout heatmap by (default: "Cell line")
#   project_name - Name for cross-sample output folder (default: "all_samples")
#
# Data paths (derived from base_path and project_name):
#   - Input:  {base_path}/scAboslute-obj/
#   - Per-sample output: {base_path}/analysis_per_sample/
#   - Cross-sample output: {base_path}/results_{project_name}/
#
# Examples:
#   ./run_pipeline.sh
#   ./run_pipeline.sh samples/PEO_samples.csv
#   ./run_pipeline.sh samples/ALLSAMPLES_metadata.csv "/Volumes/MyDrive/data" 100
#   ./run_pipeline.sh samples/ALLSAMPLES_metadata.csv "/Volumes/MyDrive/data" 100 config/qc_params_relaxed.csv
#   ./run_pipeline.sh samples/ALLSAMPLES_metadata.csv "/Volumes/MyDrive/data" 100 config/qc_params_default.csv feature1
#   ./run_pipeline.sh samples/PEO_samples.csv "/Volumes/MyDrive/data" 100 config/qc_params_default.csv "Cell line" peo
#
# ==============================================================================

SAMPLES_CSV="${1:-samples/ALLSAMPLES_metadata.csv}"
BASE_PATH="${2:-/Volumes/LenovoPS8/FIbackup/sc_analysis}"
BIN_SIZE="${3:-100}"
QC_CONFIG="${4:-config/qc_params_default.csv}"
GROUP_COL="${5:-Cell line}"
PROJECT_NAME="${6:-all_samples}"

# Derived paths
OBJ_BASE="$BASE_PATH/scAboslute-obj"
POST_BASE="$BASE_PATH/analysis_per_sample"
OUT_BASE="$BASE_PATH/results_${PROJECT_NAME}"

echo "=============================================="
echo "       Complete Analysis Pipeline"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  Samples CSV    : $SAMPLES_CSV"
echo "  Base path      : $BASE_PATH"
echo "  Bin size       : $BIN_SIZE"
echo "  QC config      : $QC_CONFIG"
echo "  Project name   : $PROJECT_NAME"
echo "  Group col      : $GROUP_COL"
echo ""
echo "Data paths:"
echo "  Input (scAbsolute objects)  : $OBJ_BASE"
echo "  Per-sample output           : $POST_BASE"
echo "  Cross-sample results        : $OUT_BASE"
echo ""
echo "=============================================="
echo ""

# Check if QC config exists
if [ ! -f "$QC_CONFIG" ]; then
    echo "ERROR: QC config file not found: $QC_CONFIG"
    echo "Available configs:"
    ls -1 config/*.csv 2>/dev/null || echo "  (none found in config/)"
    exit 1
fi

# Check if samples CSV exists
if [ ! -f "$SAMPLES_CSV" ]; then
    echo "ERROR: Samples CSV not found: $SAMPLES_CSV"
    exit 1
fi

# Check if input directory exists
if [ ! -d "$OBJ_BASE" ]; then
    echo "ERROR: Input directory not found: $OBJ_BASE"
    echo "Check your base_path argument."
    exit 1
fi

# Step 1: Generate outlier summaries
echo "Step 1: Generating outlier summaries..."
echo "----------------------------------------------"
Rscript scripts/01_generate_outlier_summaries.R \
    "$SAMPLES_CSV" \
    "$OBJ_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE" \
    "$QC_CONFIG"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed"
    exit 1
fi

# Step 2: Combine summaries
echo ""
echo "Step 2: Combining summaries..."
echo "----------------------------------------------"
Rscript scripts/02_combine_outlier_summaries.R \
    "$SAMPLES_CSV" \
    "$OUT_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed"
    exit 1
fi

# Step 3: Visualize
echo ""
echo "Step 3: Generating visualizations..."
echo "----------------------------------------------"

# Find the most recent output file from step 2
STEP2_OUTPUT=$(ls -t "$OUT_BASE"/*_*_*.csv 2>/dev/null | head -1)

if [ -z "$STEP2_OUTPUT" ]; then
    echo "ERROR: No output file found from Step 2"
    exit 1
fi

echo "Using input: $STEP2_OUTPUT"

Rscript scripts/03_visualize_summary.R \
    "$STEP2_OUTPUT" \
    "$OUT_BASE" \
    feature1 \
    Sample

if [ $? -ne 0 ]; then
    echo "ERROR: Step 3 failed"
    exit 1
fi

# Step 4: Generate dropout summaries
echo ""
echo "Step 4: Generating dropout summaries..."
echo "----------------------------------------------"
Rscript scripts/04_generate_dropout_summaries.R \
    "$SAMPLES_CSV" \
    "$OBJ_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 4 failed"
    exit 1
fi

# Step 5: Visualize dropout heatmap
echo ""
echo "Step 5: Generating dropout heatmap..."
echo "----------------------------------------------"
Rscript scripts/05_visualize_dropout_heatmap.R \
    "$SAMPLES_CSV" \
    "$POST_BASE" \
    "$OUT_BASE" \
    "$BIN_SIZE" \
    "$GROUP_COL"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 5 failed"
    exit 1
fi

# Step 6: Per-sample dropout heatmaps
echo ""
echo "Step 6: Generating per-sample dropout plots..."
echo "----------------------------------------------"
Rscript scripts/06_visualize_sample_dropouts.R \
    "$SAMPLES_CSV" \
    "$POST_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 6 failed"
    exit 1
fi

# Step 7: Generate read distribution summaries
echo ""
echo "Step 7: Generating read distribution summaries..."
echo "----------------------------------------------"
Rscript scripts/07_generate_read_distribution.R \
    "$SAMPLES_CSV" \
    "$OBJ_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 7 failed"
    exit 1
fi

echo ""
echo "=============================================="
echo "          Pipeline Complete"
echo "=============================================="
echo ""
echo "Cross-sample results : $OUT_BASE"
echo "QC summary           : $STEP2_OUTPUT"
echo "Dropout summary      : $POST_BASE/dropout_summary_combined.csv"
echo "Dropout heatmap      : $OUT_BASE/dropout_heatmap_$(date +%Y%m%d).pdf"
echo "Per-sample plots     : $POST_BASE/SLX-*/figures/"
echo "Read distribution    : $POST_BASE/read_distribution/ (summary) and $POST_BASE/SLX-*/figures/ (plots)"
echo ""
