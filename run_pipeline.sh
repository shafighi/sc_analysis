#!/bin/bash
# ==============================================================================
# Run the complete outlier summary pipeline
# ==============================================================================
#
# Usage: ./run_pipeline.sh [samples_csv] [base_path] [bin_size] [qc_config]
#
# Arguments:
#   samples_csv - CSV file with sample IDs (default: samples/ALLSAMPLES_metadata.csv)
#   base_path   - Base path to data directories (default: /Volumes/LenovoPS8/FI backup/sc_analysis)
#   bin_size    - Bin size for analysis (default: 100)
#   qc_config   - QC parameters config file (default: config/qc_params_default.csv)
#
# Data paths (derived from base_path):
#   - Input:  {base_path}/scAboslute-obj/
#   - Output: {base_path}/post-scAbsolute/
#   - Final:  {base_path}/all_samples_23july2024/output/
#
# Examples:
#   ./run_pipeline.sh
#   ./run_pipeline.sh samples/PEO_samples.csv
#   ./run_pipeline.sh samples/ALLSAMPLES_metadata.csv "/Volumes/MyDrive/data" 100
#   ./run_pipeline.sh samples/ALLSAMPLES_metadata.csv "/Volumes/MyDrive/data" 100 config/qc_params_relaxed.csv
#
# ==============================================================================

SAMPLES_CSV="${1:-samples/ALLSAMPLES_metadata.csv}"
BASE_PATH="${2:-/Volumes/LenovoPS8/FI backup/sc_analysis}"
BIN_SIZE="${3:-100}"
QC_CONFIG="${4:-config/qc_params_default.csv}"

# Derived paths
OBJ_BASE="$BASE_PATH/scAboslute-obj"
POST_BASE="$BASE_PATH/post-scAbsolute"
OUT_BASE="$BASE_PATH/all_samples_23july2024/output"

echo "=============================================="
echo "       Outlier Summary Pipeline"
echo "=============================================="
echo ""
echo "Configuration:"
echo "  Samples CSV : $SAMPLES_CSV"
echo "  Base path   : $BASE_PATH"
echo "  Bin size    : $BIN_SIZE"
echo "  QC config   : $QC_CONFIG"
echo ""
echo "Data paths:"
echo "  Input (scAbsolute objects) : $OBJ_BASE"
echo "  Intermediate (post-scAbsolute) : $POST_BASE"
echo "  Final output : $OUT_BASE"
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
STEP2_OUTPUT=$(ls -t "$OUT_BASE"/output_*.csv 2>/dev/null | head -1)

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

echo ""
echo "=============================================="
echo "          Pipeline Complete"
echo "=============================================="
echo ""
echo "Output location: $OUT_BASE"
echo "Main output file: $STEP2_OUTPUT"
echo ""
