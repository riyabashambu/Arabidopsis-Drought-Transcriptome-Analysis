#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 8: Gene Quantification using featureCounts
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

BAM_DIR="$PROJECT_DIR/results/bam"
COUNT_DIR="$PROJECT_DIR/results/counts"
GTF="$PROJECT_DIR/data/reference/Arabidopsis_thaliana.TAIR10.60.gtf"

mkdir -p "$COUNT_DIR"

echo "======================================="
echo "Running featureCounts"
echo "======================================="

featureCounts \
    -T 8 \
    -p \
    -t exon \
    -g gene_id \
    -a "$GTF" \
    -o "$COUNT_DIR/gene_counts.txt" \
    "$BAM_DIR"/*.sorted.bam

echo ""
echo "======================================="
echo "Gene counting completed!"
echo "======================================="
