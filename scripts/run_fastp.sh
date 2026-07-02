#!/bin/bash

# ==========================================================
# Arabidopsis Drought Transcriptome Analysis
# Step 3: Read Trimming using fastp
# ==========================================================

set -e

PROJECT_DIR="$HOME/Arabidopsis-Drought-Transcriptome-Analysis"

RAW_DIR="$PROJECT_DIR/data/raw_fastq"
TRIM_DIR="$PROJECT_DIR/data/trimmed_fastq"
REPORT_DIR="$PROJECT_DIR/results/fastp"

mkdir -p "$TRIM_DIR"
mkdir -p "$REPORT_DIR"

echo "======================================="
echo "Running fastp..."
echo "======================================="

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
echo "======================================="
echo "Processing $SAMPLE"
echo "======================================="

fastp \
    -i "$RAW_DIR/${SAMPLE}_1.fastq" \
    -I "$RAW_DIR/${SAMPLE}_2.fastq" \
    -o "$TRIM_DIR/${SAMPLE}_1.trimmed.fastq" \
    -O "$TRIM_DIR/${SAMPLE}_2.trimmed.fastq" \
    --detect_adapter_for_pe \
    --thread 8 \
    --qualified_quality_phred 20 \
    --length_required 30 \
    --html "$REPORT_DIR/${SAMPLE}_fastp.html" \
    --json "$REPORT_DIR/${SAMPLE}_fastp.json"

echo "$SAMPLE completed."

done

echo ""
echo "======================================="
echo "All samples trimmed successfully!"
echo "======================================="
