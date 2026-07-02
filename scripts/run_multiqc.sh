#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 2: Aggregate FastQC Reports using MultiQC
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

FASTQC_DIR="$PROJECT_DIR/results/fastqc/raw"
MULTIQC_DIR="$PROJECT_DIR/results/multiqc/raw"

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
