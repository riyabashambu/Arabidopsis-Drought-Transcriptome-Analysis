#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Download RNA-seq FASTQ files from NCBI SRA
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

RAW_DIR="$PROJECT_DIR/data/raw_fastq"
TMP_DIR="$PROJECT_DIR/data/tmp"

mkdir -p "$RAW_DIR"
mkdir -p "$TMP_DIR"

echo "======================================="
echo "Downloading RNA-seq samples..."
echo "======================================="

RUNS=(
SRR32771561
SRR32771560
SRR32771559
SRR32771552
SRR32771551
SRR32771550
)

for RUN in "${RUNS[@]}"
do
    echo ""
    echo "========================================"
    echo "Downloading $RUN"
    echo "========================================"

    prefetch "$RUN" --output-directory "$TMP_DIR"

    fasterq-dump \
        "$TMP_DIR/$RUN" \
        --split-files \
        --threads 8 \
        --temp "$TMP_DIR" \
        --outdir "$RAW_DIR"

    rm -rf "$TMP_DIR/$RUN"

    echo "$RUN finished."
done

echo ""
echo "======================================="
echo "All downloads completed!"
echo "======================================="
