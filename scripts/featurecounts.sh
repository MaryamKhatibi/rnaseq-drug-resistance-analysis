#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

BAM_DIR="$PROJECT_DIR/results/mapping"
GTF="$PROJECT_DIR/data/reference/annotation/gencode.v50.annotation.gtf"
OUTPUT_DIR="$PROJECT_DIR/results/counts"

mkdir -p "$OUTPUT_DIR"

featureCounts \
    -T 2 \
    -p \
    --countReadPairs \
    -B \
    -C \
    -s 2 \
    -a "$GTF" \
    -o "$OUTPUT_DIR/gene_counts.txt" \
    "$BAM_DIR"/*.sorted.bam

echo "featureCounts completed."
echo "Output: $OUTPUT_DIR/gene_counts.txt"