#!/bin/bash
# gatk_split_gvcf_by_chrom.sh
# Merge two cohort GVCFs (e.g. NSCP and SCAP) chromosome by chromosome:
# for every chromosome, SelectVariants extracts that chromosome from both
# GVCFs, CombineGVCFs merges them, and GenotypeGVCFs produces a VCF — all
# chromosomes processed in parallel background jobs.
#
# This avoids the memory/storage blow-up of merging whole-cohort GVCFs at once.
#
# NOTE: edit the chromosome list to match the reference genome:
#   CHM13: chr1-chr22, chrX, chrY
#   YAO  : chr1-chr22, chrX, chrY, chrM
#   hg38 : extract contig names from the FASTA header first (too many contigs
#          to hard-code) — see the companion snippet in the README.

# --- Directory layout (edit to your own paths) ------------------------------
GVCF1="/path/to/project/combine_YAO/combined_NSCP_final-YAO_raw_variants.g.vcf.gz"
GVCF2="/path/to/project/combine_YAO/combined_SCAP_final-YAO_raw_variants.g.vcf.gz"
OUTPUT_DIR="/path/to/project/combine_YAO/"
REFERENCE="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
LOG_DIR="/path/to/project/combine_YAO/logs/"
# -----------------------------------------------------------------------------

mkdir -p "$OUTPUT_DIR"
mkdir -p "$LOG_DIR"

chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12
             chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
             chrX chrY chrM)

process_chromosome() {
    chr=$1
    echo "Processing chromosome: $chr"

    # Split both GVCFs by chromosome
    gatk SelectVariants \
        -R $REFERENCE \
        -V $GVCF1 \
        -L $chr \
        -O $OUTPUT_DIR/sample1_${chr}.g.vcf.gz &> $LOG_DIR/selectvariants_sample1_${chr}.log
    if [ $? -ne 0 ]; then
        echo "Error in SelectVariants for $GVCF1, chromosome: $chr"
        cat $LOG_DIR/selectvariants_sample1_${chr}.log
        exit 1
    fi
    gatk SelectVariants \
        -R $REFERENCE \
        -V $GVCF2 \
        -L $chr \
        -O $OUTPUT_DIR/sample2_${chr}.g.vcf.gz &> $LOG_DIR/selectvariants_sample2_${chr}.log
    if [ $? -ne 0 ]; then
        echo "Error in SelectVariants for $GVCF2, chromosome: $chr"
        cat $LOG_DIR/selectvariants_sample2_${chr}.log
        exit 1
    fi

    # Combine the two GVCFs for this chromosome
    gatk CombineGVCFs \
        -R $REFERENCE \
        -V $OUTPUT_DIR/sample1_${chr}.g.vcf.gz \
        -V $OUTPUT_DIR/sample2_${chr}.g.vcf.gz \
        -O $OUTPUT_DIR/combined_${chr}.g.vcf.gz
    if [ $? -ne 0 ]; then
        echo "Error in CombineGVCFs for chromosome: $chr"
        exit 1
    fi

    # Genotype the combined GVCF into a VCF
    gatk GenotypeGVCFs \
        -R $REFERENCE \
        -V $OUTPUT_DIR/combined_${chr}.g.vcf.gz \
        -O $OUTPUT_DIR/${chr}.vcf.gz
    if [ $? -ne 0 ]; then
        echo "Error in GenotypeGVCFs for chromosome: $chr"
        exit 1
    fi

    echo "Completed processing chromosome: $chr"
}

# Export variables and function so background subshells can access them
export GVCF1
export GVCF2
export OUTPUT_DIR
export REFERENCE
export LOG_DIR
export -f process_chromosome

# Process all chromosomes in parallel, keeping track of each job's PID so
# failures inside a background subshell (exit 1) are not silently lost.
pids=()
chr_of_pid=()
for chr in "${chromosomes[@]}"; do
    process_chromosome $chr &
    pids+=($!)
    chr_of_pid+=("$chr")
done

fail=0
for i in "${!pids[@]}"; do
    if ! wait "${pids[$i]}"; then
        echo "ERROR: processing failed for chromosome ${chr_of_pid[$i]}"
        fail=1
    fi
done

if [ "$fail" -ne 0 ]; then
    echo "One or more chromosomes failed — check the logs in ${LOG_DIR}."
    exit 1
fi
echo "All tasks completed successfully."
