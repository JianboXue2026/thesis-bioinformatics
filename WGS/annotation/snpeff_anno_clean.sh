#!/bin/bash
# snpeff_anno_clean.sh
# Annotate a hard-filtered cohort VCF with SnpEff, remove warning lines that
# SnpEff prints to stdout (they corrupt the redirected VCF for tabix/bcftools),
# then compress and index the cleaned annotated VCF.
#
# Usage:
#   nohup ./snpeff_anno_clean.sh <cohort_name> <ref_tag> &
# Example:
#   nohup ./snpeff_anno_clean.sh NSCP_final YAO &

name="combined_$1"   # e.g. NSCP_final -> combined_NSCP_final

# --- Directory layout (edit to your own paths) ------------------------------
input_dir="/path/to/project/Final_vcf/03.filtered_final"
stats_dir="/path/to/project/Final_vcf/04.snpeff_anno/snpeff_report"
output_dir_raw="/path/to/project/Final_vcf/04.snpeff_anno/snpeff_raw"
output_dir_cle="/path/to/project/Final_vcf/04.snpeff_anno/snpeff_clean"
wrs_log_dir="/path/to/logs/01.Log_for_snpeff-vcf2anno-vcf"

# SnpEff database names (built locally for each reference genome)
ref_database_hg38="hg38-lo"
ref_database_CHM13="CHM13-lo"
ref_database_YAO="YAO-lo-v1.1"
# -----------------------------------------------------------------------------

input_hg38="${input_dir}/${name}-hg38_fil_var.vcf"
stats_hg38="${stats_dir}/${name}-hg38_stats-reports.html"
output_hg38_anno_raw="${output_dir_raw}/${name}-hg38_anno_raw.vcf"
output_hg38_anno_cle="${output_dir_cle}/${name}-hg38_anno_clean.vcf"

input_CHM13="${input_dir}/${name}-CHM13_fil_var.vcf"
stats_CHM13="${stats_dir}/${name}-CHM13_stats-reports.html"
output_CHM13_anno_raw="${output_dir_raw}/${name}-CHM13_anno_raw.vcf"
output_CHM13_anno_cle="${output_dir_cle}/${name}-CHM13_anno_clean.vcf"

input_YAO="${input_dir}/${name}-YAO_fil_var.vcf"
stats_YAO="${stats_dir}/${name}-YAO_stats-reports.html"
output_YAO_anno_raw="${output_dir_raw}/${name}-YAO_anno_raw.vcf"
output_YAO_anno_cle="${output_dir_cle}/${name}-YAO_anno_clean.vcf"

ref_database_var="ref_database_$2"
input_var="input_$2"
stats_var="stats_$2"
output_var_anno_raw="output_$2_anno_raw"
output_var_anno_cle="output_$2_anno_cle"
log_var="${wrs_log_dir}/log_$2_fil-vcf_2_snpeff-anno-vcf.log"

nohup bash -c "
  snpeff \
    ${!ref_database_var} \
    ${!input_var} \
    -stats ${!stats_var} \
    -o vcf \
    > ${!output_var_anno_raw}
  # Remove WARNING lines that SnpEff mixes into the redirected VCF
  awk '!/^WARNING/' ${!output_var_anno_raw} > ${!output_var_anno_cle}
  bgzip -k -@ 6 ${!output_var_anno_cle}
  tabix -p vcf ${!output_var_anno_cle}.gz
" > "${log_var}" 2>&1 &
