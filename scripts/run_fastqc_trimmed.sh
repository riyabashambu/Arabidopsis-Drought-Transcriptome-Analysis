#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 4: FastQC on Trimmed Reads
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

TRIM_DIR="$PROJECT_DIR/data/trimmed_fastq"
OUT_DIR="$PROJECT_DIR/results/fastqc/trimmed"

mkdir -p "$OUT_DIR"

echo "======================================="
echo "Running FastQC on Trimmed Reads..."
echo "======================================="

for FILE in "$TRIM_DIR"/*.trimmed.fastq
do
    echo "Processing $(basename "$FILE")"

    fastqc \
        "$FILE" \
        --threads 8 \
        --outdir "$OUT_DIR"

done

echo ""
echo "======================================="
echo "FastQC Completed!"
echo "======================================="
