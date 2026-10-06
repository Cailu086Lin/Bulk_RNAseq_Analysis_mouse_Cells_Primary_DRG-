#!/bin/bash

set -u

# Input folder containing the original FASTQ files
input_folder="/mnt/easystore/DBH_26/Stephenson/FASTQ"

# Output folder for trimmed FASTQ files
output_folder="/mnt/easystore/DBH_26/Stephenson/FASTQ/trimmed_fastq"

mkdir -p "$output_folder"

# Find all biological sample R1 FASTQs.
# This excludes Undetermined reads and avoids files already inside trimmed_fastq.
mapfile -t fastq_files < <(
    find "$input_folder" \
        -maxdepth 1 \
        -type f \
        -name "*_R1_001.fastq.gz" \
        ! -name "Undetermined*" \
        | sort
)

total=${#fastq_files[@]}

if [[ "$total" -eq 0 ]]; then
    echo "ERROR: No FASTQ files were found in:"
    echo "  $input_folder"
    exit 1
fi

echo "========================================="
echo "Starting Trim Galore on $total samples"
echo "Library type: single-end"
echo "Input:  $input_folder"
echo "Output: $output_folder"
echo "========================================="

current=0
failed=0

for R1 in "${fastq_files[@]}"; do
    current=$((current + 1))
    filename=$(basename "$R1")

    echo
    echo "[$current/$total] Processing $filename"
    echo "----------------------------------------"

    if trim_galore \
        --fastqc \
        --cores 8 \
        --output_dir "$output_folder" \
        "$R1"
    then
        echo "Successfully trimmed: $filename"
    else
        echo "ERROR trimming: $filename"
        failed=$((failed + 1))
    fi
done

echo
echo "========================================="
echo "Trim Galore processing complete"
echo "Total samples:  $total"
echo "Failed samples: $failed"
echo "Output folder:  $output_folder"
echo "========================================="

if [[ "$failed" -gt 0 ]]; then
    exit 1
fi
