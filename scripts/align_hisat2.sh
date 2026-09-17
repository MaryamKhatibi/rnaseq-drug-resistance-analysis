#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT_DIR="$PROJECT_DIR/results/trimming"
OUTPUT_DIR="$PROJECT_DIR/results/mapping"

mkdir -p "$OUTPUT_DIR"

REFERENCE="$PROJECT_DIR/data/reference/hisat2_index/grch38/genome"

SAMPLES=(
    "SRR10416712"
    "SRR10416713"
    "SRR10416714"
    "SRR10416715"
    "SRR10416716"
    "SRR10416717"
    "SRR10416718"
    "SRR10416719"
    "SRR10416720"
    "SRR10416721"
    "SRR10416722"
)

for SAMPLE in "${SAMPLES[@]}"
do
    echo "Processing ${SAMPLE}..."

    hisat2 \
        -p 2 \
        --rna-strandness RF \
        --dta \
        -x "$REFERENCE" \
        -1 "$INPUT_DIR/${SAMPLE}_R1_paired.fastq.gz" \
        -2 "$INPUT_DIR/${SAMPLE}_R2_paired.fastq.gz" \
        2> "$OUTPUT_DIR/${SAMPLE}_summary.txt" \
    | samtools sort \
        -@ 2 \
        -m 256M \
        -o "$OUTPUT_DIR/${SAMPLE}.sorted.bam"

    samtools index "$OUTPUT_DIR/${SAMPLE}.sorted.bam"

    echo "$SAMPLE processing complete"
done

echo "All samples processed."