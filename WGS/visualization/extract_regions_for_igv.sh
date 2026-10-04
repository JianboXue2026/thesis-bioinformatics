#!/bin/bash
# extract_regions_for_igv.sh
# Manual variant inspection in IGV:
#   1. define BED files around candidate variants
#      - <regions>.bam.bed : +/- 2000 bp windows (for slicing BAMs)
#      - <regions>.vcf.bed : +/- 2 bp windows (for slicing VCFs)
#   2. slice all BAMs with samtools view -L, merge per cohort, index
#   3. slice the annotated VCF with bcftools view -R
#
# BED format (tab-separated): chr <start> <end> <gene_label>

# === Step 1: BED files ========================================================
# Create e.g. regions.bam.bed:
#   chr1<TAB>778719<TAB>782719<TAB>TMEM88B
#   chr1<TAB>19207674<TAB>19208074<TAB>CRB1
#   ...
# and regions.vcf.bed:
#   chr1<TAB>780717<TAB>780721<TAB>TMEM88B
#   ...

# === Step 2: slice + merge BAMs ==============================================
cd /path/to/project/000.bam

for bam_file in *.bam; do
    output_bam="${bam_file%.bam}_extracted.bam"
    samtools view -@ 6 -h -b -L /path/to/bed/regions.bam.bed "$bam_file" \
        -o ./diff-region-bam/"$output_bam" \
        && samtools index -b ./diff-region-bam/"$output_bam"
done

cd ./diff-region-bam
# merge per cohort (sample-ID prefixes differ between cohorts)
samtools merge -@ 24 -c -p merged-NSCP.bam CAP*.bam && samtools index -b merged-NSCP.bam
samtools merge -@ 24 -c -p merged-SCAP.bam SCP*.bam && samtools index -b merged-SCAP.bam
samtools merge -@ 24 -c -p merged-all.bam merged-NSCP.bam merged-SCAP.bam \
    && samtools index -b merged-all.bam

# === Step 3: slice the annotated VCF =========================================
bcftools view -R /path/to/bed/regions.vcf.bed \
    /path/to/snpeff_anno/all-YAO-GATK-fil-snpeff-raw.vcf.gz \
    -o all-YAO-GATK-cut.vcf

# Load merged-all.bam (+ the reference genome) and the sliced VCF into IGV
# and inspect read coverage and genotypes at the candidate sites.
