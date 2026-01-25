#!/bin/bash
# Run the complete outlier summary pipeline
# Usage: ./run_pipeline.sh [samples_csv] [base_path] [bin_size]

SAMPLES_CSV="${1:-samples/ALLSAMPLES_metadata.csv}"
BASE_PATH="${2:-/Volumes/LenovoPS8/FI backup/sc_analysis}"
BIN_SIZE="${3:-100}"

OBJ_BASE="$BASE_PATH/scAboslute-obj"
POST_BASE="$BASE_PATH/post-scAbsolute"
OUT_BASE="$BASE_PATH/all_samples_23july2024/output"

echo "=== Outlier Summary Pipeline ==="
echo "Samples CSV: $SAMPLES_CSV"
echo "Base path: $BASE_PATH"
echo "Bin size: $BIN_SIZE"
echo ""

# Step 1: Generate outlier summaries
echo "Step 1: Generating outlier summaries..."
Rscript scripts/01_generate_outlier_summaries.R \
    "$SAMPLES_CSV" \
    "$OBJ_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

# Step 2: Combine summaries
echo ""
echo "Step 2: Combining summaries..."
Rscript scripts/02_combine_outlier_summaries.R \
    "$SAMPLES_CSV" \
    "$OUT_BASE" \
    "$POST_BASE" \
    "$BIN_SIZE"

# Step 3: Visualize
echo ""
echo "Step 3: Generating visualizations..."
Rscript scripts/03_visualize_summary.R \
    "$OUT_BASE/outlier_summary_table_meta_combined.csv" \
    "$OUT_BASE" \
    feature1 \
    Sample

echo ""
echo "=== Pipeline Complete ==="
echo "Output: $OUT_BASE"
