#!/bin/bash
# gatk_gvcf_to_vcf.sh
# Genotype a combined cohort GVCF into a multi-sample VCF with GATK
# GenotypeGVCFs, then compress and index the result.
#
# Usage:
#   nohup ./gatk_gvcf_to_vcf.sh <cohort_name> <ref_tag> &
# Example:
#   nohup ./gatk_gvcf_to_vcf.sh NSCP_final YAO &

echo "Start processing..."

# source activate <gatk_env>

name_1="combined_$1"   # e.g. NSCP_final -> combined_NSCP_final

# --- Directory layout (edit to your own paths) ------------------------------
ref_seq_hg38="/path/to/ref_seq/hg38/hg38.fa"
ref_seq_CHM13="/path/to/ref_seq/CHM13/T2T-CHM13v2.0.fna"
ref_seq_YAO="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
gvcf_dir="/path/to/project/GATK_gvcf"
final_vcf_dir="/path/to/project/Final_vcf"
wrs_log_dir="/path/to/logs/Log_for_gatk_gvcf"
# -----------------------------------------------------------------------------

input_hg38_1="${gvcf_dir}/${name_1}-hg38_raw_variants.g.vcf.gz"
output_hg38="${final_vcf_dir}/${name_1}-hg38_raw_variants.vcf"
input_CHM13_1="${gvcf_dir}/${name_1}-CHM13_raw_variants.g.vcf.gz"
output_CHM13="${final_vcf_dir}/${name_1}-CHM13_raw_variants.vcf"
input_YAO_1="${gvcf_dir}/${name_1}-YAO_raw_variants.g.vcf.gz"
output_YAO="${final_vcf_dir}/${name_1}-YAO_raw_variants.vcf"

ref_seq_var="ref_seq_$2"
input_1_var="input_$2_1"
output_var="output_$2"

nohup bash -c "
  gatk GenotypeGVCFs \
    -R ${!ref_seq_var} \
    -V ${!input_1_var} \
    -O ${!output_var}
  bgzip -k -@ 6 ${!output_var}
  tabix -p vcf ${!output_var}.gz
" > "${wrs_log_dir}/${name_1}_$2_gvcftovcf.log.txt" 2>&1 &

echo "GenotypeGVCFs started in the background."
