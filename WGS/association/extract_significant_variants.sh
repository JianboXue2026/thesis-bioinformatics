#!/bin/bash
# extract_significant_variants.sh
# Extract PLINK-significant variants (P < 0.05) from association results and
# pull the matching records out of an annotated VCF:
#   Step 1: .assoc -> CSV of rows with P < 0.05
#   Step 2: CSV    -> BED (chr, pos-1, pos, variant-id)
#   Step 3: BED    -> cleaned, tab-separated BED
#   Step 4: bcftools view -R extracts the variants from the annotated VCF
#
# Run inside the directory holding the .assoc files.

# === Step 1: .assoc -> P005/<prefix>.csv ====================================
mkdir -p P005
for file in *.assoc; do
    filename=$(basename "$file" .assoc)
    # column 9 = P value; keep header + significant rows
    awk 'BEGIN {FS=OFS=" "} NR==1 {print; next} $9+0 < 0.05' "$file" > "P005/${filename}.csv"
    echo "Processed $file -> P005/${filename}.csv"
done

# === Step 2: CSV -> BED ======================================================
mkdir -p BED
for file in P005/*.csv; do
    filename=$(basename "$file" .csv)
    awk 'BEGIN {FS=OFS=" "} NR>1 {print "chr"$1, $3-1, $3, $2}' "$file" > "BED/${filename}.bed"
    echo "Processed $file -> BED/${filename}.bed"
done

# === Step 3: clean a BED (tab separator, drop malformed lines) ==============
# example for one file — loop over BED/*.bed as needed
input_bed="BED/all-hg38-fil.bed"
final_bed="BED/final_cleaned_all-hg38-fil.bed"
awk 'BEGIN {FS="[ \t]+"; OFS="\t"} {print $1, $2, $3, $4}' "$input_bed" | awk 'NF == 4' > "$final_bed"
echo "Cleaned BED file -> $final_bed"

# === Step 4: extract variants from the annotated VCF ========================
# -Oz so the .vcf.gz output is actually bgzip-compressed (tabix/bcftools need it)
bcftools view -R BED/final_cleaned_all-hg38-fil.bed \
    /path/to/snpeff_anno/all-hg38-GATK-fil-snpeff-cleaned.vcf.gz \
    -Oz -o /path/to/P005.vcf/all-hg38-GATK-fil-snpeff_annoed.vcf.gz

echo "All files processed."
