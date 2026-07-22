#!/usr/bin/env bash
# Run the standalone post-scUnique visualization pipeline.

set -euo pipefail

PIPELINE_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
RSCRIPT_BIN=${RSCRIPT_BIN:-Rscript}

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: $0 <scunique_result_dir> [output_dir] [prefix]" >&2
    exit 2
fi

RESULT_DIR=$1
OUTPUT_DIR=${2:-$RESULT_DIR/post_scunique}
PREFIX=${3:-}

cmd=(
    "$RSCRIPT_BIN"
    "$PIPELINE_DIR/scripts/generate_visualizations.R"
    "$RESULT_DIR"
    "$OUTPUT_DIR"
)

if [[ -n "$PREFIX" ]]; then
    cmd+=("$PREFIX")
fi

"${cmd[@]}"
