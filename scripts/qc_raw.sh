#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

DATA_DIR="$PROJECT_DIR/data/fastq"
OUTPUT_DIR="$PROJECT_DIR/results/qc/raw"

mkdir -p "$OUTPUT_DIR"

for file in "$DATA_DIR"/*.fastq.gz
do
    echo "Running FastQC on $file"
    fastqc "$file" -o "$OUTPUT_DIR"
done

echo "FastQC analysis complete. Results are in $OUTPUT_DIR"
