#!/bin/bash
# annovar_annotate.sh
# Annotate a filtered cohort VCF with ANNOVAR (table_annovar.pl, gene-based
# refGene protocol, VCF output).
#
# The reference genome is inferred from the VCF file name
# (combined_<GROUP>_<REF>_fil_var.vcf -> <REF> = CHM13 / YAO / hg38), and the
# matching database (<REF>db) built by annovar_build_db.sh is selected.
#
# Usage:
#   screen -S annovar_anno
#   bash annovar_annotate.sh <vcf_file_name>
# Example:
#   bash annovar_annotate.sh combined_NSCP_final-CHM13_fil_var.vcf

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi
VCF_FILE=$1

# --- Directory layout (edit to your own paths) ------------------------------
VCF_DIR="/path/to/project/filtered_vcf"
DB_DIR="/path/to/project/ANNOVAR_anno/ref_database"
OUTPUT_DIR="/path/to/project/ANNOVAR_anno/annovar_raw_vcf"
ANNOVAR_DIR="/path/to/tools/annovar"
# -----------------------------------------------------------------------------

# Parse group / reference genome from the file name
FILE_NAME=$(basename "$VCF_FILE")
REF_GENOME=$(echo "$FILE_NAME" | cut -d'-' -f2 | cut -d'_' -f1)

case $REF_GENOME in
    "CHM13") DB_PATH="$DB_DIR/CHM13db" ;;
    "YAO")   DB_PATH="$DB_DIR/YAOdb"   ;;
    "hg38")  DB_PATH="$DB_DIR/hg38db"  ;;
    *)
        echo "Unknown reference genome: $REF_GENOME"
        exit 1
        ;;
esac

mkdir -p $OUTPUT_DIR
OUTPUT_FILE="$OUTPUT_DIR/${FILE_NAME%.vcf}_annovar_raw.vcf"

# Run ANNOVAR
perl $ANNOVAR_DIR/table_annovar.pl $VCF_DIR/$VCF_FILE $DB_PATH \
    -buildver $REF_GENOME \
    -out ${OUTPUT_FILE%.vcf} \
    -remove \
    -protocol refGene \
    -operation g \
    -nastring . \
    -vcfinput

if [ $? -ne 0 ]; then
    echo "Error annotating $VCF_FILE"
    exit 1
fi

# Rename the multi-annotation output to a standard VCF name
mv ${OUTPUT_FILE%.vcf}.${REF_GENOME}_multianno.vcf $OUTPUT_FILE
echo "Annotated VCF file saved as $OUTPUT_FILE"
