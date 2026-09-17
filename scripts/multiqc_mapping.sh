#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT_DIR="$PROJECT_DIR/results/mapping"
OUTPUT_DIR="$PROJECT_DIR/results/qc/mapping"

mkdir -p "$OUTPUT_DIR"

multiqc "$INPUT_DIR" -o "$OUTPUT_DIR"

echo "MultiQC mapping report generated in $OUTPUT_DIR"
