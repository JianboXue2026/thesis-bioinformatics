#!/bin/bash
# deepvariant_batch.sh
# Batch variant calling with DeepVariant (GPU Docker image, v1.6.1).
#
# The sample list is read from a CSV file (one sample per line):
#   SAMPLE_ID,REF_GENOME,SEX
#   CAP001,CHM13,1
#   CAP002,CHM13,1
#   ... (SEX: 1 = male, 0 = female; male samples get --haploid_contigs=chrX,chrY)
#
# IMPORTANT: this script pins the container to one GPU (--gpus '"device=0"').
# Running several unpinned DeepVariant containers in parallel caused GPU
# out-of-memory errors — duplicate this script with device=1 for a second GPU
# and never launch more than one container per GPU.
#
# Usage:
#   screen -S DV_GPU0
#   bash deepvariant_batch.sh samples.csv
#
# The reference FASTA must be indexed beforehand:
#   samtools faidx <ref.fa>

# Check the sample file argument
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample_file>"
    echo "Example: $0 samples.csv"
    exit 1
fi
SAMPLE_FILE=$1
if [ ! -f "$SAMPLE_FILE" ]; then
    echo "Sample file ${SAMPLE_FILE} not found, exiting..."
    exit 1
fi

# DeepVariant version
BIN_VERSION="1.6.1"

# --- Directory layout (edit to your own paths) ------------------------------
BAM_DIR="/path/to/project/000.bam"
REF_DIR="/path/to/project/VEP_ref"
VCF_OUTPUT_DIR="/path/to/project/DV_vcf"
GVCF_OUTPUT_DIR="/path/to/project/DV_gvcf"
LOG_DIR="/path/to/project/DV_log"
# -----------------------------------------------------------------------------

# Reference genome file mapping
declare -A REF_FILES
REF_FILES=(
    ["hg38"]="hg38.fa"
    ["CHM13"]="T2T-CHM13v2.0.fna"
    ["YAO"]="T2T-Yao-hp.v1.1.fasta"
)

# Read the sample file (skip the header line)
tail -n +2 "$SAMPLE_FILE" | while IFS=',' read -r SAMPLE_ID REF_GENOME SEX; do
    echo "Read from CSV: SAMPLE_ID=${SAMPLE_ID}, REF_GENOME=${REF_GENOME}, SEX=${SEX}"

    BAM_FILE="${BAM_DIR}/${SAMPLE_ID}-${REF_GENOME}.sort.fixmate.samtools_dedup.bam"
    REF_FILE="${REF_DIR}/${REF_FILES[$REF_GENOME]}"
    VCF_OUTPUT_FILE="${VCF_OUTPUT_DIR}/${SAMPLE_ID}-${REF_GENOME}.vcf"
    GVCF_OUTPUT_FILE="${GVCF_OUTPUT_DIR}/${SAMPLE_ID}-${REF_GENOME}.gvcf"

    if [ ! -f "$BAM_FILE" ]; then
        echo "BAM file for ${SAMPLE_ID} not found, skipping..."
        continue
    fi
    if [ ! -f "$REF_FILE" ]; then
        echo "Reference file for ${REF_GENOME} not found, skipping..."
        continue
    fi

    SAMPLE_LOG_DIR="${LOG_DIR}/${SAMPLE_ID}"
    mkdir -p "$SAMPLE_LOG_DIR"

    HAPLOID_CONTIGS=""
    if [ "$SEX" -eq 1 ]; then
        HAPLOID_CONTIGS="--haploid_contigs=chrX,chrY"
    fi

    echo "Processing sample ${SAMPLE_ID} with reference genome ${REF_GENOME}..."
    docker run --gpus '"device=0"' \
      -v "${BAM_DIR}":"/input" \
      -v "${VCF_OUTPUT_DIR}:/output" \
      -v "${GVCF_OUTPUT_DIR}:/output_gvcf" \
      -v "${LOG_DIR}:/output/logs" \
      -v "${REF_DIR}:/input/ref" \
      google/deepvariant:"${BIN_VERSION}-gpu" \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type=WGS \
      --ref=/input/ref/${REF_FILES[$REF_GENOME]} \
      --reads=/input/${SAMPLE_ID}-${REF_GENOME}.sort.fixmate.samtools_dedup.bam \
      --output_vcf=/output/${SAMPLE_ID}-${REF_GENOME}.vcf \
      --output_gvcf=/output_gvcf/${SAMPLE_ID}-${REF_GENOME}.gvcf \
      --num_shards=48 \
      --logging_dir=/output/logs/${SAMPLE_ID} \
      $HAPLOID_CONTIGS \
      --dry_run=false > "${SAMPLE_LOG_DIR}/deepvariant.log" 2>&1
done
