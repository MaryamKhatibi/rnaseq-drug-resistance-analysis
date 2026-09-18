#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="$PROJECT_DIR/results/counts/gene_counts.txt"
OUTPUT="$PROJECT_DIR/results/counts/count_matrix.txt"

mkdir -p "$(dirname "$OUTPUT")"

awk 'BEGIN {OFS="\t"}
NR == 1 {next}
NR == 2 {
    print $1, \
          "BxPC3_S1", "BxPC3_S2", "BxPC3_S3", \
          "BxPC3_GR1", "BxPC3_GR2", "BxPC3_GR3", \
          "CFPAC1_S1", "CFPAC1_S2", "CFPAC1_S3", \
          "CFPAC1_GR1", "CFPAC1_GR2", "CFPAC1_GR3"
    next
}
{
    print $1, $7, $8, $9, $10, $11, $12, \
          $13, $14, $15, $16, $17, $18
}' "$INPUT" > "$OUTPUT"

echo "Count matrix prepared."
echo "Input : $INPUT"
echo "Output: $OUTPUT"