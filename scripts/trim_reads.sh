#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT_DIR="$PROJECT_DIR/data/fastq"
OUTPUT_DIR="$PROJECT_DIR/results/trimming"

mkdir -p "$OUTPUT_DIR"

ADAPTERS="$CONDA_PREFIX/share/trimmomatic-0.41-0/adapters/TruSeq3-PE.fa"

for R1 in "$INPUT_DIR"/*_1.fastq.gz
do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="$INPUT_DIR/${SAMPLE}_2.fastq.gz"

    echo "Trimming $SAMPLE"

    trimmomatic PE -threads 4 \
        "$R1" \
        "$R2" \
        "$OUTPUT_DIR/${SAMPLE}_R1_paired.fastq.gz" \
        "$OUTPUT_DIR/${SAMPLE}_R1_unpaired.fastq.gz" \
        "$OUTPUT_DIR/${SAMPLE}_R2_paired.fastq.gz" \
        "$OUTPUT_DIR/${SAMPLE}_R2_unpaired.fastq.gz" \
        ILLUMINACLIP:"$ADAPTERS":2:30:10 \
        MINLEN:36

done

echo "Trimming complete. Results are in $OUTPUT_DIR"
