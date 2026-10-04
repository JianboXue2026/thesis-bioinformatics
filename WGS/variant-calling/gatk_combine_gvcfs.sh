#!/bin/bash
# gatk_combine_gvcfs.sh
# Hierarchically merge per-sample GVCFs into cohort-level combined GVCFs with
# GATK CombineGVCFs, then compress and index the result.
#
# Merging is performed level by level (5 samples -> lv1 parts, lv1 parts ->
# lv2 parts, ... -> <cohort>_final) because whole-cohort GVCFs are extremely
# large. The reference genome is selected via the last argument (hg38 / CHM13
# / YAO) using indirect variable expansion.
#
# Usage:
#   nohup ./gatk_combine_gvcfs.sh <in1> <in2> <in3> <in4> <out_name> <ref_tag> &
# Example (input names are GVCF name prefixes, e.g. "lv2_part01-05"):
#   nohup ./gatk_combine_gvcfs.sh lv2_part01-05 lv2_part06-10 lv2_part11-15 lv1_part16 NSCP_final YAO &

echo "Start processing..."
core_processing=4

# source activate <gatk_env>

name_1="combined_$1"
name_2="combined_$2"
name_3="combined_$3"
name_4="combined_$4"
output_name=$5

# --- Directory layout (edit to your own paths) ------------------------------
ref_seq_hg38="/path/to/ref_seq/hg38/hg38.fa"
ref_seq_CHM13="/path/to/ref_seq/CHM13/T2T-CHM13v2.0.fna"
ref_seq_YAO="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
gvcf_dir="/path/to/project/GATK_gvcf"
wrs_log_dir="/path/to/logs/Log_for_gatk_gvcf"
# -----------------------------------------------------------------------------

input_hg38_1="${gvcf_dir}/${name_1}-hg38_raw_variants.g.vcf.gz"
input_hg38_2="${gvcf_dir}/${name_2}-hg38_raw_variants.g.vcf.gz"
input_hg38_3="${gvcf_dir}/${name_3}-hg38_raw_variants.g.vcf.gz"
input_hg38_4="${gvcf_dir}/${name_4}-hg38_raw_variants.g.vcf.gz"
output_hg38="${gvcf_dir}/combined_${output_name}-hg38_raw_variants.g.vcf"

input_CHM13_1="${gvcf_dir}/${name_1}-CHM13_raw_variants.g.vcf.gz"
input_CHM13_2="${gvcf_dir}/${name_2}-CHM13_raw_variants.g.vcf.gz"
input_CHM13_3="${gvcf_dir}/${name_3}-CHM13_raw_variants.g.vcf.gz"
input_CHM13_4="${gvcf_dir}/${name_4}-CHM13_raw_variants.g.vcf.gz"
output_CHM13="${gvcf_dir}/combined_${output_name}-CHM13_raw_variants.g.vcf"

input_YAO_1="${gvcf_dir}/${name_1}-YAO_raw_variants.g.vcf.gz"
input_YAO_2="${gvcf_dir}/${name_2}-YAO_raw_variants.g.vcf.gz"
input_YAO_3="${gvcf_dir}/${name_3}-YAO_raw_variants.g.vcf.gz"
input_YAO_4="${gvcf_dir}/${name_4}-YAO_raw_variants.g.vcf.gz"
output_YAO="${gvcf_dir}/combined_${output_name}-YAO_raw_variants.g.vcf"

# Select the reference genome via indirect expansion: ref_seq_<tag>, input_<tag>_N ...
ref_seq_var="ref_seq_$6"
input_1_var="input_$6_1"
input_2_var="input_$6_2"
input_3_var="input_$6_3"
input_4_var="input_$6_4"
output_var="output_$6"

nohup bash -c "
  gatk CombineGVCFs \
    -R ${!ref_seq_var} \
    -V ${!input_1_var} \
    -V ${!input_2_var} \
    -V ${!input_3_var} \
    -V ${!input_4_var} \
    -O ${!output_var}
  bgzip -k -@ ${core_processing} ${!output_var}
  tabix -p vcf ${!output_var}.gz
" > "${wrs_log_dir}/combined_${output_name}_$6_gvcfgz.log" 2>&1 &

echo "CombineGVCFs started in the background."
