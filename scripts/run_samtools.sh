#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 7: Process SAM files using SAMtools
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

ALIGN_DIR="$PROJECT_DIR/results/alignment"
BAM_DIR="$PROJECT_DIR/results/bam"

mkdir -p "$BAM_DIR"

SAMPLES=(
SRR32771561
SRR32771560
SRR32771559
SRR32771552
SRR32771551
SRR32771550
)

echo "======================================="
echo "Processing alignment files..."
echo "======================================="

for SAMPLE in "${SAMPLES[@]}"
do

    echo ""
    echo "Processing $SAMPLE"

    # Convert SAM → BAM
    samtools view \
        -@ 8 \
        -bS "$ALIGN_DIR/${SAMPLE}.sam" \
        > "$BAM_DIR/${SAMPLE}.bam"

    # Sort BAM
    samtools sort \
        -@ 8 \
        "$BAM_DIR/${SAMPLE}.bam" \
        -o "$BAM_DIR/${SAMPLE}.sorted.bam"

    # Index BAM
    samtools index \
        "$BAM_DIR/${SAMPLE}.sorted.bam"

    # Alignment statistics
    samtools flagstat \
        "$BAM_DIR/${SAMPLE}.sorted.bam" \
        > "$BAM_DIR/${SAMPLE}_flagstat.txt"

    # Remove unsorted BAM to save space
    rm "$BAM_DIR/${SAMPLE}.bam"

    echo "$SAMPLE completed."

done

echo ""
echo "======================================="
echo "SAMtools processing completed!"
echo "======================================="
