#!/bin/bash
# workflow_fastq_to_dedup_bam.sh
# Fully automated workflow from clean FASTQ to deduplicated, indexed,
# statistics-reported BAM — one call per sample, three reference genomes in
# parallel background jobs.
#
# This variant also supports sample renaming (old_name -> new_name), used to
# standardise sample IDs (e.g. raw sequencing IDs -> unified study IDs).
#
# Usage:
#   nohup ./workflow_fastq_to_dedup_bam.sh <old_name> <new_name> &
# For samples that need no renaming, pass the same ID twice:
#   nohup ./workflow_fastq_to_dedup_bam.sh CAP001 CAP001 &

echo "Start processing..."
core_processing=8
core_converting=8

old_name=$1
new_name=$2

# --- Directory layout (edit to your own paths) -----------------------------
cleaned_data_dir="/path/to/project/Cleaned_Data"
ref_seq_dir="/path/to/ref_seq"
wrs_temp_dir="/path/to/temp"
bma_hg38_dir="/path/to/project/BMA/HG38"
bma_chm13_dir="/path/to/project/BMA/CHM13"
bma_yao_dir="/path/to/project/BMA/T2T-YAO"
samtools_stats_dir="/path/to/project/Samtools_stats"
wrs_log_dir="/path/to/logs/Log_for_workflow"
# ----------------------------------------------------------------------------

input_1="${cleaned_data_dir}/${old_name}_1.fq.gz"
input_2="${cleaned_data_dir}/${old_name}_2.fq.gz"
rg_string="@RG\tID:${new_name}\tSM:${new_name}\tLB:WGS\tPL:ILLUMINA"

# --- hg38 -------------------------------------------------------------------
nohup bash -c "
  bwa mem -t ${core_processing} -R '${rg_string}' ${ref_seq_dir}/hg38/hg38.fa ${input_1} ${input_2} | \
    samtools sort -n -m 3G -@ ${core_converting} -T ${wrs_temp_dir}/${new_name}-hg38.temp.bwa -o ${bma_hg38_dir}/${new_name}-hg38.bam -
  samtools fixmate -@ ${core_processing} -m ${bma_hg38_dir}/${new_name}-hg38.bam ${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.bam
  samtools sort -@ ${core_processing} -T ${wrs_temp_dir}/${new_name}-hg38.temp.dedup -m 3G ${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.bam | \
    samtools markdup -@ ${core_processing} - ${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.samtools_dedup.bam
  samtools index -@ ${core_processing} ${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.samtools_dedup.bam
  samtools stats -@ ${core_processing} ${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.samtools_dedup.bam > ${samtools_stats_dir}/${new_name}-hg38-stat.txt
  # Uncomment to delete intermediates after verification:
  # rm \"${bma_hg38_dir}/${new_name}-hg38.bam\" \"${bma_hg38_dir}/${new_name}-hg38.sort.fixmate.bam\"
" > "${wrs_log_dir}/${new_name}_hg38.log" 2>&1 &

# --- CHM13 ------------------------------------------------------------------
nohup bash -c "
  bwa mem -t ${core_processing} -R '${rg_string}' ${ref_seq_dir}/CHM13/T2T-CHM13v2.0.fna ${input_1} ${input_2} | \
    samtools sort -n -m 3G -@ ${core_converting} -T ${wrs_temp_dir}/${new_name}-CHM13.temp.bwa -o ${bma_chm13_dir}/${new_name}-CHM13.bam -
  samtools fixmate -@ ${core_processing} -m ${bma_chm13_dir}/${new_name}-CHM13.bam ${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.bam
  samtools sort -@ ${core_processing} -T ${wrs_temp_dir}/${new_name}-CHM13.temp.dedup -m 3G ${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.bam | \
    samtools markdup -@ ${core_processing} - ${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.samtools_dedup.bam
  samtools index -@ ${core_processing} ${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.samtools_dedup.bam
  samtools stats -@ ${core_processing} ${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.samtools_dedup.bam > ${samtools_stats_dir}/${new_name}-CHM13-stat.txt
  # rm \"${bma_chm13_dir}/${new_name}-CHM13.bam\" \"${bma_chm13_dir}/${new_name}-CHM13.sort.fixmate.bam\"
" > "${wrs_log_dir}/${new_name}_CHM13.log" 2>&1 &

# --- T2T-YAO ----------------------------------------------------------------
nohup bash -c "
  bwa mem -t ${core_processing} -R '${rg_string}' ${ref_seq_dir}/YAO/T2T-Yao-hp.v1.1.fasta ${input_1} ${input_2} | \
    samtools sort -n -m 3G -@ ${core_converting} -T ${wrs_temp_dir}/${new_name}-YAO.temp.bwa -o ${bma_yao_dir}/${new_name}-YAO.bam -
  samtools fixmate -@ ${core_processing} -m ${bma_yao_dir}/${new_name}-YAO.bam ${bma_yao_dir}/${new_name}-YAO.sort.fixmate.bam
  samtools sort -@ ${core_processing} -T ${wrs_temp_dir}/${new_name}-YAO.temp.dedup -m 3G ${bma_yao_dir}/${new_name}-YAO.sort.fixmate.bam | \
    samtools markdup -@ ${core_processing} - ${bma_yao_dir}/${new_name}-YAO.sort.fixmate.samtools_dedup.bam
  samtools index -@ ${core_processing} ${bma_yao_dir}/${new_name}-YAO.sort.fixmate.samtools_dedup.bam
  samtools stats -@ ${core_processing} ${bma_yao_dir}/${new_name}-YAO.sort.fixmate.samtools_dedup.bam > ${samtools_stats_dir}/${new_name}-YAO-stat.txt
  # rm \"${bma_yao_dir}/${new_name}-YAO.bam\" \"${bma_yao_dir}/${new_name}-YAO.sort.fixmate.bam\"
" > "${wrs_log_dir}/${new_name}_YAO.log" 2>&1 &

# Wait for all background jobs
wait
echo "All processes finished."
