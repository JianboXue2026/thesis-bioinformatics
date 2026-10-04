#!/bin/bash
# gatk_haplotypecaller_gvcf.sh
# Per-sample variant calling with GATK HaplotypeCaller in GVCF mode, followed
# by bgzip compression and tabix indexing, for the three reference genomes.
#
# Prerequisites per reference genome (run once):
#   gzip -d hg38.fa.gz            # GATK cannot read .gz reference directly
#   samtools faidx <ref.fa>       # -> <ref.fa>.fai
#   samtools dict <ref.fa> > <ref.dict>
#
# Usage (run inside screen or with nohup so it survives SSH disconnects):
#   nohup ./gatk_haplotypecaller_gvcf.sh <sample_id> &

echo "Start processing..."
core_processing=4
name=$1

# Activate an environment providing GATK 4.x, e.g.:
# source activate <gatk_env>

# --- Directory layout (edit to your own paths) ------------------------------
ref_seq_hg38="/path/to/ref_seq/hg38/hg38.fa"
ref_seq_CHM13="/path/to/ref_seq/CHM13/T2T-CHM13v2.0.fna"
ref_seq_YAO="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
gvcf_dir="/path/to/project/GATK_gvcf"
bam_hg38="/path/to/project/BMA/HG38/${name}-hg38.sort.fixmate.samtools_dedup.bam"
bam_CHM13="/path/to/project/BMA/CHM13/${name}-CHM13.sort.fixmate.samtools_dedup.bam"
bam_YAO="/path/to/project/BMA/T2T-YAO/${name}-YAO.sort.fixmate.samtools_dedup.bam"
wrs_log_dir="/path/to/logs/Log_for_gatk_gvcf"
# -----------------------------------------------------------------------------

# --- hg38 ---
nohup bash -c "
  gatk HaplotypeCaller \
    -R ${ref_seq_hg38} \
    -I ${bam_hg38} \
    -ERC GVCF \
    -O ${gvcf_dir}/${name}-hg38_raw_variants.g.vcf
  bgzip -k -@ ${core_processing} ${gvcf_dir}/${name}-hg38_raw_variants.g.vcf
  tabix -p vcf ${gvcf_dir}/${name}-hg38_raw_variants.g.vcf.gz
" > "${wrs_log_dir}/${name}_hg38_gvcfgz.log" 2>&1 &

# --- CHM13 ---
nohup bash -c "
  gatk HaplotypeCaller \
    -R ${ref_seq_CHM13} \
    -I ${bam_CHM13} \
    -ERC GVCF \
    -O ${gvcf_dir}/${name}-CHM13_raw_variants.g.vcf
  bgzip -k -@ ${core_processing} ${gvcf_dir}/${name}-CHM13_raw_variants.g.vcf
  tabix -p vcf ${gvcf_dir}/${name}-CHM13_raw_variants.g.vcf.gz
" > "${wrs_log_dir}/${name}_CHM13_gvcfgz.log" 2>&1 &

# --- YAO ---
nohup bash -c "
  gatk HaplotypeCaller \
    -R ${ref_seq_YAO} \
    -I ${bam_YAO} \
    -ERC GVCF \
    -O ${gvcf_dir}/${name}-YAO_raw_variants.g.vcf
  bgzip -k -@ ${core_processing} ${gvcf_dir}/${name}-YAO_raw_variants.g.vcf
  tabix -p vcf ${gvcf_dir}/${name}-YAO_raw_variants.g.vcf.gz
" > "${wrs_log_dir}/${name}_YAO_gvcfgz.log" 2>&1 &

wait
echo "All processes finished."
