#!/bin/bash
# gatk_merge_chr_vcfs.sh
# Merge the per-chromosome VCFs produced by gatk_split_gvcf_by_chrom.sh into
# one cohort-level VCF with GATK MergeVcfs.

# --- Directory layout (edit to your own paths) ------------------------------
input_dir="/path/to/project/combine_YAO"
output_file="/path/to/project/combine_YAO/combined-all-YAO.vcf.gz"
ref_genome="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
# -----------------------------------------------------------------------------

# Collect all per-chromosome VCF inputs
input_files=$(ls $input_dir/chr*.vcf.gz)

# Build the -I parameter list
input_params=""
for file in $input_files;
do
    input_params+=" -I $file"
done

# Merge
gatk MergeVcfs $input_params -O $output_file -R $ref_genome

# Verify the output file exists
if [ -f $output_file ]; then
    echo "VCF files merged successfully into $output_file"
else
    echo "Error: merging VCF files failed"
fi
