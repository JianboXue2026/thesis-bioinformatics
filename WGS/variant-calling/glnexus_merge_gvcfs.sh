#!/bin/bash
# glnexus_merge_gvcfs.sh
# Joint-call DeepVariant GVCFs with GLnexus and convert the merged BCF to a
# compressed VCF. GLnexus is equivalent to `gatk CombineGVCFs + GenotypeGVCFs`
# but much faster for DeepVariant output.
#
# Two configurations are used in this study:
#   DeepVariantWGS        — merge + standard filtering (default DeepVariant WGS model)
#   DeepVariant_unfiltered — merge WITHOUT filtering (keeps raw variant calls)
#
# NOTE: GLnexus writes its ".GLnexus.DB" temp directory into the current
# working directory. Run each glnexus_cli instance in its OWN temp directory —
# the temp files of concurrent runs collide if they share a directory.
#
# Input file list format: one .gvcf path per line, e.g.
#   ./CAP001-CHM13.gvcf
#   ./CAP002-CHM13.gvcf

# --- User-configurable section ----------------------------------------------
combine_name="NSCP-CHM13"                                  # output name prefix
config="DeepVariantWGS"                                    # or DeepVariant_unfiltered
list_file="/path/to/project/DV_gvcf/file_list/${combine_name}.txt"
output_bcf="/path/to/project/DV_gvcf/${combine_name}.bcf"
output_vcf_gz="/path/to/project/DV_gvcf/${combine_name}.vcf.gz"
temp_dir="/path/to/project/DV_gvcf/temp/${combine_name}"  # dedicated temp dir
# -----------------------------------------------------------------------------

mkdir -p "$temp_dir"
cd "$temp_dir"

# Merge the GVCFs
glnexus_cli --config ${config} \
            --mem-gbytes 100 \
            --threads 18 \
            --list "${list_file}" \
            > "${output_bcf}"

# Convert BCF to bgzipped VCF
bcftools view "${output_bcf}" | bgzip -@ 6 -c > "${output_vcf_gz}"
