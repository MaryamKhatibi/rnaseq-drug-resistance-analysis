#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT_DIR="$PROJECT_DIR/results/trimming"
OUTPUT_DIR="$PROJECT_DIR/results/qc/trimmed"

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*_paired.fastq.gz
do
    echo "Running FastQC on $file"
    fastqc "$file" -o "$OUTPUT_DIR"
done

echo "FastQC analysis complete. Results are in $OUTPUT_DIR"
