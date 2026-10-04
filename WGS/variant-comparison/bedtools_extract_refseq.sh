#!/bin/bash
# bedtools_extract_refseq.sh
# Extract reference sequences around candidate variant sites with bedtools
# getfasta — used for primer / motif inspection and downstream validation.
#
# BED format (tab-separated): chr <start> <end> [gene_label]
#   +/- 500 bp windows for focused sequence extraction
#   +/- 200 bp windows for gene-level sequence extraction

# --- User-configurable section ------------------------------------------------
ref_fasta="/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
bed_500="/path/to/bed/regions_500bp.bed"     # e.g. chr1  197827817  197828817
bed_200="/path/to/bed/regions_200bp.bed"     # e.g. chr1  192078541  192078941  CRB1
out_500="regions_500bp.fasta"
out_200="regions_200bp_genes.fasta"
# -----------------------------------------------------------------------------

bedtools getfasta -fi ${ref_fasta} -bed ${bed_500} -fo ${out_500}
bedtools getfasta -fi ${ref_fasta} -bed ${bed_200} -fo ${out_200}

echo "Extracted sequences: ${out_500}, ${out_200}"
