#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 5: MultiQC on Trimmed Reads
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

FASTQC_DIR="$PROJECT_DIR/results/fastqc/trimmed"
MULTIQC_DIR="$PROJECT_DIR/results/multiqc/trimmed"

mkdir -p "$MULTIQC_DIR"

echo "======================================="
echo "Running MultiQC..."
echo "======================================="

multiqc \
    "$FASTQC_DIR" \
    --outdir "$MULTIQC_DIR"

echo ""
echo "======================================="
echo "MultiQC Completed!"
echo "======================================="
