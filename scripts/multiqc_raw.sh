#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT_DIR="$PROJECT_DIR/results/qc/raw"
OUTPUT_DIR="$PROJECT_DIR/results/qc/multiqc"

mkdir -p "$OUTPUT_DIR"

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC analysis complete. Report is in $OUTPUT_DIR"

