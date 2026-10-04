#!/bin/bash
# plink_qc.sh
# Quality control of a PLINK bed file set:
#   --mind 0.05   remove samples with > 5% missing genotypes
#   --geno 0.05   remove variants with > 5% missing genotypes
#   --maf  0.01   remove variants with MAF < 1%
#   --hwe  1e-6   remove variants significantly deviating from HWE (p < 1e-6)
#
# NOTE: for very-low-frequency variant sets the thresholds may need relaxing
# (e.g. --mind 0.60 --geno 0.60) — pick per dataset and document the choice.
#
# Usage:
#   screen -S PLINK_fil
#   bash plink_qc.sh <input_prefix>

# Check the input prefix argument
if [ -z "$1" ]; then
  echo "Usage: $0 <input_prefix>"
  exit 1
fi

# --- Directory layout (edit to your own paths) ------------------------------
input_prefix=$1
input_dir="/path/to/GATK-PLINK/01.trans2bed"
output_dir="/path/to/GATK-PLINK/02.filtered_bed"
# -----------------------------------------------------------------------------

output_prefix="${output_dir}/${input_prefix}-fil"
mkdir -p "${output_dir}"

plink2 --bfile "${input_dir}/${input_prefix}" \
       --mind 0.05 \
       --geno 0.05 \
       --maf 0.01 \
       --hwe 1e-6 \
       --make-bed \
       --allow-extra-chr \
       --out "${output_prefix}"

echo "Quality control completed. Output files are saved in ${output_dir}"
