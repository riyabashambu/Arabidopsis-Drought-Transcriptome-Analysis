#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 1: Quality Control with FastQC
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

RAW_DIR="$PROJECT_DIR/data/raw_fastq"
OUT_DIR="$PROJECT_DIR/results/fastqc/raw"

mkdir -p "$OUT_DIR"

echo "======================================="
echo "Running FastQC..."
echo "======================================="

for FILE in "$RAW_DIR"/*.fastq
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
echo "Reports saved in:"
echo "$OUT_DIR"
echo "======================================="
