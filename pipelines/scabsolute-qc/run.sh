#!/usr/bin/env bash
# Run the complete post-scAbsolute QC and reporting pipeline.

set -euo pipefail

PIPELINE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
CALL_DIR=$PWD
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}

usage() {
    cat <<'EOF'
Usage:
  run.sh <samples_csv> <base_path> [bin_size] [qc_config] [group_col] [project_name]

Arguments:
  samples_csv   CSV containing a required 'sample' column
  base_path     Working data directory
  bin_size      Bin size used in scAbsolute filenames (default: 100)
  qc_config     QC parameter CSV (default: config/qc_params_default.csv)
  group_col     Metadata column for grouped plots (default: Cell line)
  project_name  Cross-sample output suffix (default: all_samples)

Expected input:
  <base_path>/scAboslute-obj/SLX-<sample>_<bin_size>.rds
EOF
}

if [[ $# -lt 2 || $# -gt 6 ]]; then
    usage >&2
    exit 2
fi

absolute_from_call_dir() {
    case "$1" in
        /*) printf '%s\n' "$1" ;;
        *)  printf '%s/%s\n' "$CALL_DIR" "$1" ;;
    esac
}

SAMPLES_CSV=$(absolute_from_call_dir "$1")
BASE_PATH=$(absolute_from_call_dir "$2")
BIN_SIZE=${3:-100}
QC_CONFIG=${4:-$PIPELINE_DIR/config/qc_params_default.csv}
GROUP_COL=${5:-Cell line}
PROJECT_NAME=${6:-all_samples}

if [[ "$QC_CONFIG" != /* ]]; then
    QC_CONFIG=$(absolute_from_call_dir "$QC_CONFIG")
fi

OBJ_BASE=$BASE_PATH/scAboslute-obj
POST_BASE=$BASE_PATH/analysis_per_sample
OUT_BASE=$BASE_PATH/results_$PROJECT_NAME

if [[ ! -f "$SAMPLES_CSV" ]]; then
    echo "ERROR: sample manifest not found: $SAMPLES_CSV" >&2
    exit 1
fi
if [[ ! -f "$QC_CONFIG" ]]; then
    echo "ERROR: QC config not found: $QC_CONFIG" >&2
    exit 1
fi
if [[ ! -d "$OBJ_BASE" ]]; then
    echo "ERROR: scAbsolute input directory not found: $OBJ_BASE" >&2
    exit 1
fi

mkdir -p "$POST_BASE" "$OUT_BASE"
cd "$PIPELINE_DIR"

echo "Post-scAbsolute QC pipeline"
echo "  samples : $SAMPLES_CSV"
echo "  inputs  : $OBJ_BASE"
echo "  per-cell: $POST_BASE"
echo "  summary : $OUT_BASE"
echo "  bin size: $BIN_SIZE"
echo "  config  : $QC_CONFIG"

echo
echo "[1/7] Per-sample QC classification and profile reports"
"$RSCRIPT_BIN" scripts/01_generate_outlier_summaries.R \
    "$SAMPLES_CSV" "$OBJ_BASE" "$POST_BASE" "$BIN_SIZE" "$QC_CONFIG"

echo
echo "[2/7] Combined QC summary"
"$RSCRIPT_BIN" scripts/02_combine_outlier_summaries.R \
    "$SAMPLES_CSV" "$OUT_BASE" "$POST_BASE" "$BIN_SIZE"

STEP2_OUTPUT=""
for candidate in "$OUT_BASE"/*.csv; do
    [[ -e "$candidate" ]] || continue
    if [[ -z "$STEP2_OUTPUT" || "$candidate" -nt "$STEP2_OUTPUT" ]]; then
        STEP2_OUTPUT=$candidate
    fi
done
if [[ -z "$STEP2_OUTPUT" ]]; then
    echo "ERROR: step 2 did not produce a summary CSV in $OUT_BASE" >&2
    exit 1
fi

echo
echo "[3/7] Cross-sample QC visualizations"
"$RSCRIPT_BIN" scripts/03_visualize_summary.R \
    "$STEP2_OUTPUT" "$OUT_BASE" "$GROUP_COL" Sample

echo
echo "[4/7] Per-sample dropout summaries"
"$RSCRIPT_BIN" scripts/04_generate_dropout_summaries.R \
    "$SAMPLES_CSV" "$OBJ_BASE" "$POST_BASE" "$BIN_SIZE"

echo
echo "[5/7] Cross-sample dropout heatmap"
"$RSCRIPT_BIN" scripts/05_visualize_dropout_heatmap.R \
    "$SAMPLES_CSV" "$POST_BASE" "$OUT_BASE" "$BIN_SIZE" "$GROUP_COL"

echo
echo "[6/7] Per-sample dropout plots"
"$RSCRIPT_BIN" scripts/06_visualize_sample_dropouts.R \
    "$SAMPLES_CSV" "$POST_BASE" "$POST_BASE" "$BIN_SIZE"

echo
echo "[7/7] Read-distribution summaries"
"$RSCRIPT_BIN" scripts/07_generate_read_distribution.R \
    "$SAMPLES_CSV" "$OBJ_BASE" "$POST_BASE" "$BIN_SIZE"

echo
echo "Pipeline complete"
echo "  QC summary        : $STEP2_OUTPUT"
echo "  Cross-sample plots: $OUT_BASE"
echo "  Per-sample outputs: $POST_BASE"
