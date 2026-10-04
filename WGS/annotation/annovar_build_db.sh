#!/bin/bash
# annovar_build_db.sh
# Build a custom ANNOVAR gene-definition database from a reference genome's
# GFF3 + FASTA:
#   1) gff3ToGenePred  : GFF3 -> GenePred format (refGene.txt)
#   2) retrieve_seq_from_fasta.pl : transcript FASTA (refGeneMrna.fa)
#
# Prerequisite: gff3ToGenePred binary (UCSC) available on PATH.
# For CHM13, run fix_gff3_order.sh on the GFF3 first; for YAO, use an
# AGAT-cleaned GFF3 (see fix_gff3_parent_refs.py for the raw-GFF3 Parent issues).
#
# Usage:
#   screen -S annovar_db
#   bash annovar_build_db.sh <gff3_file> <fasta_file> <output_subdir>
# Example:
#   bash annovar_build_db.sh hg38.gff3 hg38.fa hg38db

# Check arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <gff3_file> <fasta_file> <output_subdir>"
    exit 1
fi
GFF3_FILE=$1
FASTA_FILE=$2
OUTPUT_SUBDIR=$3

# --- Directory layout (edit to your own paths) ------------------------------
REFSEQ_DIR="/path/to/project/VEP_ref"
OUTPUT_DIR="/path/to/project/ANNOVAR_anno/ref_database/$OUTPUT_SUBDIR"
ANNOVAR_SCRIPT_DIR="/path/to/tools/annovar"
# -----------------------------------------------------------------------------

mkdir -p $OUTPUT_DIR

# 1) GFF3 -> GenePred
echo "Converting GFF3 to GenePred..."
gff3ToGenePred $REFSEQ_DIR/$GFF3_FILE $OUTPUT_DIR/refGene.txt
if [ $? -ne 0 ]; then
    echo "Error: failed to convert GFF3 to GenePred."
    exit 1
fi

# 2) Generate transcript FASTA for ANNOVAR
echo "Generating mRNA FASTA file for ANNOVAR..."
perl $ANNOVAR_SCRIPT_DIR/retrieve_seq_from_fasta.pl \
    --format refGene \
    --seqfile $REFSEQ_DIR/$FASTA_FILE \
    $OUTPUT_DIR/refGene.txt \
    --out $OUTPUT_DIR/refGeneMrna.fa
if [ $? -ne 0 ]; then
    echo "Error: failed to generate mRNA FASTA file."
    exit 1
fi

echo "Done. GenePred and mRNA FASTA files are located in $OUTPUT_DIR."
echo "Rename them to <ref>_refGene.txt / <ref>_refGeneMrna.fa before annotation."
