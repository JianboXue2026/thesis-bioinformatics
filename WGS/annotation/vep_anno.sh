#!/bin/bash
# vep_anno.sh
# Annotate a hard-filtered cohort VCF with Ensembl VEP (Docker image), using
# local FASTA + GFF3 files so no internet connection or cache download is
# required.
#
# Reference FASTA/GFF3 preparation (once per genome):
#   - rename files so <name>.fasta / <name>.sorted.gff3 match the script
#   - bgzip -k <ref>.fasta ; samtools faidx <ref>.fasta.gz
#   - bgzip -k <ref>.sorted.gff3 ; tabix -p gff <ref>.sorted.gff3.gz
#   (GFF3 must be sorted by chrom, start beforehand — see fix_gff3_order.sh)
#
# NOTE: with --fasta and --gff specified, VEP runs in offline mode
# automatically. Do NOT add --offline — it triggers a cache check that fails
# with "No cache found for homo_sapiens".
#
# Usage:
#   screen -S vep_anno
#   bash vep_anno.sh <reference_genome_name> <vcf_file_name>
# Example:
#   bash vep_anno.sh T2T-Yao-hp.v1.1 combined_NSCP_final-YAO_fil_var.vcf

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reference_genome_name> <vcf_file_name>"
    exit 1
fi
REFERENCE_GENOME_NAME=$1
VCF_FILE_NAME=$2

# --- Directory layout (edit to your own paths) ------------------------------
REF_DIR="/path/to/project/VEP_ref"
VCF_DIR="/path/to/project/filtered_vcf"
OUTPUT_DIR="/path/to/project/VEP_annoed_vcf"
TEMP_OUTPUT_DIR="/opt/vep/temp_output"
# -----------------------------------------------------------------------------

REFERENCE_GENOME_PATH="${REF_DIR}/${REFERENCE_GENOME_NAME}.fasta"
GFF_FILE_PATH="${REF_DIR}/${REFERENCE_GENOME_NAME}.sorted.gff3.gz"
VCF_FILE_PATH="${VCF_DIR}/${VCF_FILE_NAME}"

# Check inputs exist
if [ ! -f "${REFERENCE_GENOME_PATH}.gz" ]; then
    echo "Reference genome file not found: ${REFERENCE_GENOME_PATH}.gz"
    exit 1
fi
if [ ! -f "${GFF_FILE_PATH}" ]; then
    echo "GFF file not found: ${GFF_FILE_PATH}"
    exit 1
fi
if [ ! -f "${VCF_FILE_PATH}" ]; then
    echo "VCF file not found: ${VCF_FILE_PATH}"
    exit 1
fi

OUTPUT_FILE_NAME="${VCF_FILE_NAME%_fil_var.vcf}_vep_anno_raw.vcf"
mkdir -p ${OUTPUT_DIR}

# --user keeps the container at the host user's permission level so the
# output directory stays writable by the host.
docker run -t -i --user $(id -u):$(id -g) \
    -v ${REF_DIR}:/opt/vep/.vep \
    -v ${VCF_DIR}:/opt/vep/input \
    -v ${OUTPUT_DIR}:${TEMP_OUTPUT_DIR} \
    ensemblorg/ensembl-vep \
    vep --input_file /opt/vep/input/${VCF_FILE_NAME} \
        --output_file ${TEMP_OUTPUT_DIR}/${OUTPUT_FILE_NAME} \
        --fasta /opt/vep/.vep/${REFERENCE_GENOME_NAME}.fasta.gz \
        --gff /opt/vep/.vep/${REFERENCE_GENOME_NAME}.sorted.gff3.gz \
        --species homo_sapiens \
        --force_overwrite \
        --vcf \
        --dir_cache /opt/vep/.vep \
        --dir_plugins /opt/vep/.vep

echo "Annotation completed. Output saved to ${OUTPUT_DIR}/${OUTPUT_FILE_NAME}"
