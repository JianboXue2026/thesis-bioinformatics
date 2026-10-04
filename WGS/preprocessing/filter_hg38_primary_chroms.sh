#!/bin/bash
# filter_hg38_primary_chroms.sh
# The UCSC hg38 FASTA contains many unplaced / alt contigs. This script keeps
# only the primary chromosomes (chr1-chr22, chrX, chrY, chrM) from the hg38
# reference FASTA and its GFF3 annotation.
#
# Usage: run inside the directory holding hg38.fa and the (sorted) GFF3.

# Chromosomes to keep
chromosomes=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9"
             "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17"
             "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" "chrM")

# --- FASTA: extract chromosome sequences into one file ---
> hg38.filtered.fa
for chr in "${chromosomes[@]}"; do
    samtools faidx hg38.fa "$chr" >> hg38.filtered.fa
done

# --- GFF3: keep records of the same chromosomes ---
> hg38.filtered.gff3
for chr in "${chromosomes[@]}"; do
    awk -v chr="$chr" '$1 == chr' hg38.sorted.gff3 >> hg38.filtered.gff3
done

# --- GFF3: sort by chromosome then start position ---
sort -k1,1 -k4,4n hg38.filtered.gff3 -o hg38.filtered.sorted.gff3

# Results: hg38.filtered.fa / hg38.filtered.gff3 / hg38.filtered.sorted.gff3
