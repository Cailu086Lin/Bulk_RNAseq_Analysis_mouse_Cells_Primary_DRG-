#!/usr/bin/env bash

set -euo pipefail

###############################################################################
# GRCm39 decoy-aware Salmon quantification
#
# Samples:
#   A1–A3 = ESP treatment
#   B1–B3 = SGE treatment
#   C1–C3 = Control
#
# Library:
#   Single-end, approximately 100-bp reads
###############################################################################

# ---------------------------------------------------------------------------
# 1. Software
# ---------------------------------------------------------------------------

SALMON="salmon"

# Confirm Salmon is available
command -v "$SALMON" >/dev/null 2>&1 || {
    echo "ERROR: Salmon was not found in PATH."
    exit 1
}

"$SALMON" --version


# ---------------------------------------------------------------------------
# 2. Directories
# ---------------------------------------------------------------------------

PROJECT_DIR="/mnt/easystore/DBH_26/Stephenson"

REF_DIR="${PROJECT_DIR}/reference/GRCm39"

FASTQ_DIR="${PROJECT_DIR}/FASTQ/trimmed_fastq"

# Change this if your trimmed FASTQs are in a different directory:
# FASTQ_DIR="${PROJECT_DIR}/trimmed_fastq"

INDEX_DIR="${REF_DIR}/salmon_index_decoy"

QUANT_DIR="${PROJECT_DIR}/salmon_GRCm39_decoy_quant"

mkdir -p "$REF_DIR"
mkdir -p "$INDEX_DIR"
mkdir -p "$QUANT_DIR"


# ---------------------------------------------------------------------------
# 3. Reference files
# ---------------------------------------------------------------------------

GENOME_GZ="${REF_DIR}/GRCm39.primary_assembly.genome.fa.gz"
TRANSCRIPTS_GZ="${REF_DIR}/gencode.vM39.transcripts.fa.gz"

DECOYS="${REF_DIR}/decoys.txt"
GENTROME_GZ="${REF_DIR}/gencode.vM39.GRCm39.gentrome.fa.gz"


# ---------------------------------------------------------------------------
# 4. Verify reference files
# ---------------------------------------------------------------------------

if [[ ! -s "$GENOME_GZ" ]]; then
    echo "ERROR: Genome FASTA was not found:"
    echo "  $GENOME_GZ"
    exit 1
fi

if [[ ! -s "$TRANSCRIPTS_GZ" ]]; then
    echo "ERROR: Transcript FASTA was not found:"
    echo "  $TRANSCRIPTS_GZ"
    exit 1
fi

echo "Genome:"
echo "  $GENOME_GZ"

echo "Transcripts:"
echo "  $TRANSCRIPTS_GZ"


# ---------------------------------------------------------------------------
# 5. Generate the chromosome/contig decoy list
# ---------------------------------------------------------------------------

echo
echo "Generating decoy list..."

gzip -cd "$GENOME_GZ" |
    grep '^>' |
    sed 's/^>//' |
    cut -d ' ' -f 1 \
    > "$DECOYS"

if [[ ! -s "$DECOYS" ]]; then
    echo "ERROR: Failed to generate decoys.txt"
    exit 1
fi

echo "Number of genome decoy sequences:"
wc -l "$DECOYS"

echo "First decoy names:"
head "$DECOYS"


# ---------------------------------------------------------------------------
# 6. Create the gentrome FASTA
#
# Transcript sequences must come first, followed by the genome.
# ---------------------------------------------------------------------------

echo
echo "Creating transcriptome + genome gentrome FASTA..."

{
    gzip -cd "$TRANSCRIPTS_GZ"
    gzip -cd "$GENOME_GZ"
} | gzip -c > "$GENTROME_GZ"

if [[ ! -s "$GENTROME_GZ" ]]; then
    echo "ERROR: Failed to create the gentrome FASTA."
    exit 1
fi

echo "Gentrome FASTA created:"
echo "  $GENTROME_GZ"


# ---------------------------------------------------------------------------
# 7. Build the GRCm39 decoy-aware Salmon index
# ---------------------------------------------------------------------------

echo
echo "Building GRCm39 decoy-aware Salmon index..."

# Remove an incomplete previous index, if present
if [[ -d "$INDEX_DIR" ]]; then
    rm -rf "$INDEX_DIR"
fi

mkdir -p "$INDEX_DIR"

"$SALMON" index \
    --transcripts "$GENTROME_GZ" \
    --decoys "$DECOYS" \
    --index "$INDEX_DIR" \
    --kmerLen 31 \
    --threads 12

if [[ ! -s "${INDEX_DIR}/versionInfo.json" ]]; then
    echo "ERROR: Salmon index generation did not complete properly."
    exit 1
fi

echo
echo "Index successfully built:"
echo "  $INDEX_DIR"


# ---------------------------------------------------------------------------
# 8. Locate and quantify all samples
# ---------------------------------------------------------------------------

samples=(
    A1 A2 A3
    B1 B2 B3
    C1 C2 C3
)

echo
echo "FASTQ directory:"
echo "  $FASTQ_DIR"

if [[ ! -d "$FASTQ_DIR" ]]; then
    echo "ERROR: FASTQ directory does not exist:"
    echo "  $FASTQ_DIR"
    exit 1
fi


for sample in "${samples[@]}"; do

    echo
    echo "================================================================="
    echo "Processing sample: $sample"
    echo "================================================================="

    # Search common Trim Galore single-end filenames.
    # The -maxdepth 1 restriction avoids accidentally finding output files
    # in nested directories.
    mapfile -t matches < <(
        find "$FASTQ_DIR" -maxdepth 1 -type f \
        \( \
            -name "${sample}_trimmed.fq.gz" -o \
            -name "${sample}_trimmed.fastq.gz" -o \
            -name "${sample}"'*_trimmed.fq.gz' -o \
            -name "${sample}"'*_trimmed.fastq.gz' -o \
            -name "${sample}"'*R1*.fq.gz' -o \
            -name "${sample}"'*R1*.fastq.gz' \
        \) | sort
    )

    if [[ ${#matches[@]} -eq 0 ]]; then
        echo "ERROR: No trimmed FASTQ found for sample $sample"
        echo "Files currently present in:"
        echo "  $FASTQ_DIR"
        ls -lh "$FASTQ_DIR"
        exit 1
    fi

    if [[ ${#matches[@]} -gt 1 ]]; then
        echo "ERROR: More than one possible FASTQ was found for $sample:"
        printf '  %s\n' "${matches[@]}"
        echo "Please make the filename pattern more specific."
        exit 1
    fi

    FASTQ="${matches[0]}"
    SAMPLE_OUT="${QUANT_DIR}/${sample}_quant"

    echo "Input FASTQ:"
    echo "  $FASTQ"

    echo "Output directory:"
    echo "  $SAMPLE_OUT"

    rm -rf "$SAMPLE_OUT"

    "$SALMON" quant \
        --index "$INDEX_DIR" \
        --libType A \
        --unmatedReads "$FASTQ" \
        --threads 8 \
        --validateMappings \
        --seqBias \
        --fldMean 200 \
        --fldSD 80 \
        --gcSizeSamp 0 \
        --output "$SAMPLE_OUT"

    if [[ ! -s "${SAMPLE_OUT}/quant.sf" ]]; then
        echo "ERROR: quant.sf was not created for $sample"
        exit 1
    fi

done


# ---------------------------------------------------------------------------
# 9. Create a mapping-rate summary
# ---------------------------------------------------------------------------

SUMMARY="${QUANT_DIR}/GRCm39_decoy_mapping_summary.tsv"

printf "sample\tprocessed_reads\tmapped_reads\tmapping_rate_percent\n" \
    > "$SUMMARY"

for sample in "${samples[@]}"; do

    META="${QUANT_DIR}/${sample}_quant/aux_info/meta_info.json"

    if command -v jq >/dev/null 2>&1; then

        processed=$(jq -r '.num_processed' "$META")
        mapped=$(jq -r '.num_mapped' "$META")
        rate=$(awk -v m="$mapped" -v n="$processed" \
            'BEGIN {printf "%.2f", 100*m/n}')

    else

        processed=$(grep -o '"num_processed":[[:space:]]*[0-9]*' "$META" |
            grep -o '[0-9]*')

        mapped=$(grep -o '"num_mapped":[[:space:]]*[0-9]*' "$META" |
            grep -o '[0-9]*')

        rate=$(awk -v m="$mapped" -v n="$processed" \
            'BEGIN {printf "%.2f", 100*m/n}')

    fi

    printf "%s\t%s\t%s\t%s\n" \
        "$sample" "$processed" "$mapped" "$rate" \
        >> "$SUMMARY"

done


echo
echo "================================================================="
echo "All GRCm39 decoy-aware Salmon quantifications completed."
echo "================================================================="

echo
echo "Quantification directory:"
echo "  $QUANT_DIR"

echo
echo "Mapping summary:"
column -t -s $'\t' "$SUMMARY" 2>/dev/null || cat "$SUMMARY"

