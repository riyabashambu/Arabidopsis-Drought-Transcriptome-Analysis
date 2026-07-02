#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 6: HISAT2 Alignment
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

INDEX="$PROJECT_DIR/data/reference/index/TAIR10_index"

TRIM_DIR="$PROJECT_DIR/data/trimmed_fastq"

ALIGN_DIR="$PROJECT_DIR/results/alignment"

mkdir -p "$ALIGN_DIR"

echo "========================================"
echo "Starting HISAT2 Alignment"
echo "========================================"

SAMPLES=(
SRR32771561
SRR32771560
SRR32771559
SRR32771552
SRR32771551
SRR32771550
)

for SAMPLE in "${SAMPLES[@]}"
do

    echo ""
    echo "========================================"
    echo "Aligning $SAMPLE"
    echo "========================================"

    hisat2 \
        -p 8 \
        -x "$INDEX" \
        -1 "$TRIM_DIR/${SAMPLE}_1.trimmed.fastq" \
        -2 "$TRIM_DIR/${SAMPLE}_2.trimmed.fastq" \
        --summary-file "$ALIGN_DIR/${SAMPLE}_summary.txt" \
        -S "$ALIGN_DIR/${SAMPLE}.sam"

    echo "$SAMPLE completed."

done

echo ""
echo "========================================"
echo "All alignments completed!"
echo "========================================"
