#!/bin/bash
# vcf_to_plink_bed.sh
# Prepare cohort VCFs for PLINK association analysis:
#   1. keep biallelic SNPs only (bcftools view -m2 -M2 -v snps)
#   2. merge the NSCP and SCAP cohort VCFs per reference genome (bcftools merge)
#   3. convert VCF -> PLINK bed/bim/fam with plink2 (--make-bed)
#
# Notes:
#   - PLINK cannot handle multiallelic sites — hence step 1.
#   - --psam supplies sex information (1 = male, 2 = female, 0 = unknown);
#     generate the .psam with csv_to_psam.py.
#   - --split-par sets the chrX PAR boundaries (values from PLINK 2 docs;
#     T2T-YAO has no published boundary so the CHM13 values were reused).
#   - hg38 carries extra contigs: add --allow-extra-chr.

# === Step 1: biallelic SNP filtering ========================================
# per-cohort VCFs
bcftools view -m2 -M2 -v snps /path/to/filtered_vcf/combined_NSCP_final-CHM13_fil_var.vcf.gz \
    -Oz -o /path/to/GATK-PLINK/00.filtered_single/NSCP-CHM13_fil.vcf.gz
bcftools view -m2 -M2 -v snps /path/to/filtered_vcf/combined_SCAP_final-CHM13_fil_var.vcf.gz \
    -Oz -o /path/to/GATK-PLINK/00.filtered_single/SCAP-CHM13_fil.vcf.gz
# (repeat for hg38 / YAO; the same applies to the DeepVariant-derived VCFs)

# === Step 2: merge NSCP + SCAP per reference genome =========================
input_dir="/path/to/project/filtered_vcf"
output_dir="/path/to/GATK-PLINK/00.filtered_single"

# bgzip + index all filtered VCFs
for vcf in "$input_dir"/*.vcf; do
    bgzip -c "$vcf" > "$vcf.gz"
    tabix -p vcf "$vcf.gz"
done

bcftools merge -Oz -o merged-CHM13.vcf.gz $(ls "$input_dir"/*-CHM13.vcf.gz)
bcftools merge -Oz -o merged-hg38.vcf.gz  $(ls "$input_dir"/*-hg38.vcf.gz)
bcftools merge -Oz -o merged-YAO.vcf.gz   $(ls "$input_dir"/*-YAO.vcf.gz)

# keep biallelic SNPs only in the merged files
bcftools view -m2 -M2 -v snps merged-CHM13.vcf.gz -Oz -o "$output_dir/merged-CHM13-single.vcf.gz"
bcftools view -m2 -M2 -v snps merged-hg38.vcf.gz  -Oz -o "$output_dir/merged-hg38-single.vcf.gz"
bcftools view -m2 -M2 -v snps merged-YAO.vcf.gz   -Oz -o "$output_dir/merged-YAO-single.vcf.gz"

# === Step 3: VCF -> PLINK bed ===============================================
# CHM13 (and YAO, reusing CHM13 PAR boundaries)
plink2 --vcf /path/to/GATK-PLINK/00.filtered_single/NSCP-CHM13_fil.vcf.gz \
    --psam /path/to/GATK-PLINK/01.trans2bed/NSCP_sex_info.psam \
    --split-par 2394410 153925835 \
    --make-bed --out /path/to/GATK-PLINK/01.trans2bed/NSCP-CHM13-gatk2plink
plink2 --vcf /path/to/GATK-PLINK/00.filtered_single/merged-CHM13-single.vcf.gz \
    --psam /path/to/GATK-PLINK/01.trans2bed/merged_sex_info.psam \
    --split-par 2394410 153925835 \
    --make-bed --out /path/to/GATK-PLINK/01.trans2bed/merged-CHM13
plink2 --vcf /path/to/GATK-PLINK/00.filtered_single/merged-YAO-single.vcf.gz \
    --psam /path/to/GATK-PLINK/01.trans2bed/merged_sex_info.psam \
    --split-par 2394410 153925835 \
    --make-bed --out /path/to/GATK-PLINK/01.trans2bed/merged-YAO

# hg38 (extra contigs -> --allow-extra-chr)
plink2 --vcf /path/to/GATK-PLINK/00.filtered_single/merged-hg38-single.vcf.gz \
    --psam /path/to/GATK-PLINK/01.trans2bed/merged_sex_info.psam \
    --split-par 2781479 155701383 \
    --make-bed --allow-extra-chr \
    --out /path/to/GATK-PLINK/01.trans2bed/merged-hg38

# DeepVariant-derived VCFs: add --vcf-half-call m to handle half-calls
plink2 --vcf /path/to/DV-PLINK/00.filtered_single/all-CHM13-DV_single.vcf.gz \
    --psam /path/to/GATK-PLINK/01.trans2bed/merged_sex_info.psam \
    --vcf-half-call m --split-par 2394410 153925835 \
    --make-bed --out /path/to/DV-PLINK/01.trans2bed/all-CHM13-DV
