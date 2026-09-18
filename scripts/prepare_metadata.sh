#!/usr/bin/env bash

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

INPUT="$PROJECT_DIR/data/metadata/SraRunTable.csv"
OUTPUT="$PROJECT_DIR/data/metadata/sample_metadata.tsv"

python - "$INPUT" "$OUTPUT" <<'PY'
import csv
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

replicate_counter = {
    "BxPC-3_Sensitive": 0,
    "BxPC-3-GR_Resistant": 0,
    "CFPAC-1_Sensitive": 0,
    "CFPAC-1-GR_Resistant": 0,
}

with open(input_file, newline="", encoding="utf-8") as infile:
    reader = csv.DictReader(infile)

    with open(output_file, "w", newline="", encoding="utf-8") as outfile:
        writer = csv.writer(outfile, delimiter="\t")

        writer.writerow([
            "sample_id",
            "geo_accession",
            "sra_run",
            "cell_line",
            "resistance",
            "replicate"
        ])

        for row in reader:
            run = row["Run"]
            geo = row["GEO_Accession (exp)"]
            cell = row["cell_line"]
            resistance = row["drug_sensitivity"]

            key = f"{cell}_{resistance}"

            if key not in replicate_counter:
                raise ValueError(
                    f"Unexpected sample: {run}\t{cell}\t{resistance}"
                )

            replicate_counter[key] += 1
            replicate = replicate_counter[key]

            if cell == "BxPC-3" and resistance == "Sensitive":
                sample_id = f"BxPC3_S{replicate}"
            elif cell == "BxPC-3-GR" and resistance == "Resistant":
                sample_id = f"BxPC3_GR{replicate}"
            elif cell == "CFPAC-1" and resistance == "Sensitive":
                sample_id = f"CFPAC1_S{replicate}"
            elif cell == "CFPAC-1-GR" and resistance == "Resistant":
                sample_id = f"CFPAC1_GR{replicate}"
            else:
                raise ValueError(
                    f"Unexpected sample: {run}\t{cell}\t{resistance}"
                )

            writer.writerow([
                sample_id,
                geo,
                run,
                cell,
                resistance,
                replicate
            ])

print("Metadata prepared.")
print(f"Input : {input_file}")
print(f"Output: {output_file}")
PY