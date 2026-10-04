#!/bin/bash
# gatk_hardfilter_vcf.sh
# Hard-filter a cohort VCF:
#   1) split into SNP and INDEL VCFs (SelectVariants)
#   2) mark failing variants with GATK-recommended hard-filter thresholds
#      (VariantFiltration) — separate parameter sets for SNP and INDEL
#   3) keep only variants that PASS (SelectVariants --exclude-filtered)
#   4) merge the filtered SNP + INDEL VCFs back together (MergeVcfs)
#
# Filter rationale (GATK best practices):
#   - single-expression filters are avoided so sites missing an annotation
#     (e.g. MQRankSum/ReadPosRankSum at non-heterozygous sites) are not
#     dropped by compound expressions.
#
# Usage:
#   nohup ./gatk_hardfilter_vcf.sh <cohort_name> <ref_tag> &
# Example:
#   nohup ./gatk_hardfilter_vcf.sh NSCP_final YAO &

echo "Start processing..."
name="combined_$1"   # e.g. NSCP_final -> combined_NSCP_final

# source activate <gatk_env>

# --- Directory layout (edit to your own paths) ------------------------------
ref_seq_hg38="/path/to/ref_seq/hg38/hg38.fa"
ref_seq_CHM13="/path/to/ref_seq/CHM13/T2T-CHM13v2.0.fna"
ref_seq_YAO="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
input_dir="/path/to/project/Final_vcf"
output_sep_raw_dir="/path/to/project/Final_vcf/01.raw_sep_snp-indel"
output_sep_fil_dir="/path/to/project/Final_vcf/02.filtered_sep_snp-indel"
output_com_fil_dir="/path/to/project/Final_vcf/03.filtered_final"
wrs_log_dir="/path/to/logs/00.Log_for_vcf2fil-vcf"
# -----------------------------------------------------------------------------

input_hg38="${input_dir}/${name}-hg38_raw_variants.vcf.gz"
output_hg38_sep_raw_snp="${output_sep_raw_dir}/${name}-hg38_raw_var-snp.vcf"
output_hg38_sep_raw_ind="${output_sep_raw_dir}/${name}-hg38_raw_var-ind.vcf"
output_hg38_sep_fil_snp_marked="${output_sep_fil_dir}/${name}-hg38_fil_var-snp-marked.vcf"
output_hg38_sep_fil_snp="${output_sep_fil_dir}/${name}-hg38_fil_var-snp.vcf"
output_hg38_sep_fil_ind_marked="${output_sep_fil_dir}/${name}-hg38_fil_var-ind-marked.vcf"
output_hg38_sep_fil_ind="${output_sep_fil_dir}/${name}-hg38_fil_var-ind.vcf"
output_hg38_com_fil="${output_com_fil_dir}/${name}-hg38_fil_var.vcf"

input_CHM13="${input_dir}/${name}-CHM13_raw_variants.vcf.gz"
output_CHM13_sep_raw_snp="${output_sep_raw_dir}/${name}-CHM13_raw_var-snp.vcf"
output_CHM13_sep_raw_ind="${output_sep_raw_dir}/${name}-CHM13_raw_var-ind.vcf"
output_CHM13_sep_fil_snp_marked="${output_sep_fil_dir}/${name}-CHM13_fil_var-snp-marked.vcf"
output_CHM13_sep_fil_snp="${output_sep_fil_dir}/${name}-CHM13_fil_var-snp.vcf"
output_CHM13_sep_fil_ind_marked="${output_sep_fil_dir}/${name}-CHM13_fil_var-ind-marked.vcf"
output_CHM13_sep_fil_ind="${output_sep_fil_dir}/${name}-CHM13_fil_var-ind.vcf"
output_CHM13_com_fil="${output_com_fil_dir}/${name}-CHM13_fil_var.vcf"

input_YAO="${input_dir}/${name}-YAO_raw_variants.vcf.gz"
output_YAO_sep_raw_snp="${output_sep_raw_dir}/${name}-YAO_raw_var-snp.vcf"
output_YAO_sep_raw_ind="${output_sep_raw_dir}/${name}-YAO_raw_var-ind.vcf"
output_YAO_sep_fil_snp_marked="${output_sep_fil_dir}/${name}-YAO_fil_var-snp-marked.vcf"
output_YAO_sep_fil_snp="${output_sep_fil_dir}/${name}-YAO_fil_var-snp.vcf"
output_YAO_sep_fil_ind_marked="${output_sep_fil_dir}/${name}-YAO_fil_var-ind-marked.vcf"
output_YAO_sep_fil_ind="${output_sep_fil_dir}/${name}-YAO_fil_var-ind.vcf"
output_YAO_com_fil="${output_com_fil_dir}/${name}-YAO_fil_var.vcf"

ref_seq_var="ref_seq_$2"
input_var="input_$2"
output_var_sep_raw_snp="output_$2_sep_raw_snp"
output_var_sep_raw_ind="output_$2_sep_raw_ind"
output_var_sep_fil_snp_marked="output_$2_sep_fil_snp_marked"
output_var_sep_fil_snp="output_$2_sep_fil_snp"
output_var_sep_fil_ind_marked="output_$2_sep_fil_ind_marked"
output_var_sep_fil_ind="output_$2_sep_fil_ind"
output_var_com_fil_dir="output_$2_com_fil"
log_var="${wrs_log_dir}/log_$2_raw-vcf_2_fil-vcf.log"

nohup bash -c "
  # 1) split raw VCF into raw SNP VCF and raw INDEL VCF
  gatk SelectVariants \
    -R ${!ref_seq_var} \
    -V ${!input_var} \
    --select-type-to-include SNP \
    -O ${!output_var_sep_raw_snp}
  gatk SelectVariants \
    -R ${!ref_seq_var} \
    -V ${!input_var} \
    --select-type-to-include INDEL \
    -O ${!output_var_sep_raw_ind}

  # 2) hard-filter the SNP VCF
  gatk VariantFiltration \
    -V ${!output_var_sep_raw_snp} \
    -filter 'QD < 2.0' --filter-name 'QD2' \
    -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
    -filter 'SOR > 3.0' --filter-name 'SOR3' \
    -filter 'FS > 60.0' --filter-name 'FS60' \
    -filter 'MQ < 40.0' --filter-name 'MQ40' \
    -filter 'MQRankSum < -12.5' --filter-name 'MQRankSum-12.5' \
    -filter 'ReadPosRankSum < -8.0' --filter-name 'ReadPosRankSum-8' \
    -O ${!output_var_sep_fil_snp_marked}
  # 3) keep PASS variants only
  gatk SelectVariants \
    -V ${!output_var_sep_fil_snp_marked} \
    --exclude-filtered \
    -O ${!output_var_sep_fil_snp}

  # 2) hard-filter the INDEL VCF
  gatk VariantFiltration \
    -V ${!output_var_sep_raw_ind} \
    -filter 'QD < 2.0' --filter-name 'QD2' \
    -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
    -filter 'FS > 200.0' --filter-name 'FS200' \
    -filter 'ReadPosRankSum < -20.0' --filter-name 'ReadPosRankSum-20' \
    -O ${!output_var_sep_fil_ind_marked}
  # 3) keep PASS variants only
  gatk SelectVariants \
    -V ${!output_var_sep_fil_ind_marked} \
    --exclude-filtered \
    -O ${!output_var_sep_fil_ind}

  # 4) merge filtered SNP + INDEL VCFs into the final filtered VCF
  gatk MergeVcfs \
    -I ${!output_var_sep_fil_snp} \
    -I ${!output_var_sep_fil_ind} \
    -O ${!output_var_com_fil_dir}
" > "${log_var}" 2>&1 &

echo "Hard filtering started in the background."
