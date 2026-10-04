#!/bin/bash
# plink_assoc.sh
# PLINK 1.9 case/control association analysis for the SCAP cohort study.
#
# Combinations run in the study (per reference genome and per caller):
#   - QC-filtered (-fil)   vs unfiltered (-nofil) bed files
#   - --assoc    (no covariates)
#   - --logistic (covariates allowed)
#   - --adjust   adds multiple-testing corrected P values to the output
#
# Key output files:
#   *.assoc / *.assoc.logistic               — per-variant statistics
#   *.assoc.adjusted / *.assoc.logistic.adjusted — corrected P values
#     UNADJ (raw), GC (genomic control), BONF, HOLM, SIDAK_SS, SIDAK_SD,
#     FDR_BH (Benjamini-Hochberg), FDR_BY (Benjamini-Yekutieli)
#   The columns of interest are typically UNADJ, GC and FDR_BH.

# --- Directory layout (edit to your own paths) ------------------------------
gatk_dir="/path/to/GATK-PLINK"
dv_dir="/path/to/DV-PLINK"
pheno_all="/path/to/GATK-PLINK/04.PLINK1-all/PLINK_phenotype.txt"
# -----------------------------------------------------------------------------

# === GATK-derived bed files ==================================================
# QC-filtered + --assoc + --adjust
for ref in CHM13 hg38 YAO; do
    plink --bfile "${gatk_dir}/02.filtered_bed/merged-${ref}-fil" \
        --pheno "${pheno_all}" --assoc --adjust --allow-extra-chr \
        --out "${gatk_dir}/03.GWAS/merged-${ref}-fil"
    # QC-filtered + --logistic + --adjust
    plink --bfile "${gatk_dir}/02.filtered_bed/merged-${ref}-fil" \
        --pheno "${pheno_all}" --logistic --adjust --allow-extra-chr \
        --out "${gatk_dir}/03.GWAS/merged-${ref}-fil-logistic"
    # unfiltered + --assoc + --adjust
    plink --bfile "${gatk_dir}/01.trans2bed/merged-${ref}" \
        --pheno "${pheno_all}" --assoc --adjust --allow-extra-chr \
        --out "${gatk_dir}/03.GWAS/merged-${ref}-nofil"
done

# === all-sample bed files (GATK) =============================================
for ref in CHM13 hg38 YAO; do
    plink --bfile "${gatk_dir}/02.filtered_bed/all-${ref}-fil" \
        --pheno "${pheno_all}" --assoc --adjust --allow-extra-chr \
        --out "${gatk_dir}/04.PLINK1-all/all-${ref}-fil"
    plink --bfile "${gatk_dir}/02.filtered_bed/all-${ref}-fil" \
        --pheno "${pheno_all}" --logistic --adjust --allow-extra-chr \
        --out "${gatk_dir}/04.PLINK1-all/all-${ref}-fil-logistic"
done

# === DeepVariant-derived bed files ===========================================
for ref in CHM13 hg38 YAO; do
    plink --bfile "${dv_dir}/02.filtered_bed/all-${ref}-DV-fil" \
        --pheno "${pheno_all}" --assoc --adjust --allow-extra-chr \
        --out "${dv_dir}/04.PLINK1-all/all-${ref}-DV-fil"
    plink --bfile "${dv_dir}/01.trans2bed/all-${ref}-DV" \
        --pheno "${pheno_all}" --assoc --adjust --allow-extra-chr \
        --out "${dv_dir}/04.PLINK1-all/all-${ref}-DV-nofil"
done
