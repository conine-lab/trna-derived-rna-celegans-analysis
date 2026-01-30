#!/usr/bin/env bash
set -euo pipefail

# Filter BAM files by reference name, index them, and compute per-base coverage
#
# Usage:
#   ./filter_bams_and_coverage.sh <input_bam_dir> <output_dir> <reference_pattern>
#
# Example:
#   ./filter_bams_and_coverage.sh \
#       tRNA/sense \
#       tRNA/sense/filtered/GluCTC \
#       GluCTC
#
# Requirements:
#   samtools >= 1.10
#   deeptools (bamCoverage)

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_bam_dir> <output_dir> <reference_pattern>"
    exit 1
fi

bam_dir="$1"
output_dir="$2"
pattern="$3"

mkdir -p "$output_dir"

for bam_file in "$bam_dir"/*.bam; do
    file_name=$(basename "$bam_file" .bam)

    filtered_bam="$output_dir/${file_name}_filtered.bam"
    filtered_bai="$output_dir/${file_name}_filtered.bai"
    coverage_out="$output_dir/${file_name}_summary_unnorm.tsv"

    echo "[filter] $bam_file (pattern: $pattern)"

    # Filter alignments by reference name (reference column match)
    samtools view -h "$bam_file" \
        | awk -v pat="$pattern" 'BEGIN{OFS="\t"} /^@/ || $3 ~ pat {print}' \
        | samtools view -b - > "$filtered_bam"

    # Index filtered BAM
    samtools index "$filtered_bam" "$filtered_bai"

    # Compute per-base coverage
    bamCoverage \
        -b "$filtered_bam" \
        -o "$coverage_out" \
        --binSize 1 \
        --outFileFormat bedgraph

done

echo "[done]"
