#!/bin/bash
# fix_gff3_order.sh
# Fix a GFF3 file whose records are not sorted (tabix refuses to index it with
# "Unsorted positions" errors): split off the header, sort the body by
# chromosome and start position, and reassemble.
#
# Used for the T2T-CHM13 annotation GFF3 before bgzip + tabix indexing and
# before ANNOVAR database building.

# --- User-configurable section ----------------------------------------------
GFF3_FILE="/path/to/project/VEP_ref/T2T-CHM13v2.0.gff3"
OUTPUT_FILE="/path/to/project/VEP_ref/T2T-CHM13v2.0.sorted.gff3"
# -----------------------------------------------------------------------------

# Temp files
HEADER_FILE=$(mktemp)
BODY_FILE=$(mktemp)

# Separate header lines and body
awk 'BEGIN {header=1}
     /^#/ {if (header) print > "'$HEADER_FILE'"}
     !/^#/ {header=0; print > "'$BODY_FILE'"}' $GFF3_FILE

# Sort the body by chromosome, then start position.
# -k1,1V gives natural contig order (chr1, chr2, ... chr10, ...) so the file
# matches the reference sequence order expected by tabix.
sort -k1,1V -k4,4n $BODY_FILE > ${BODY_FILE}.sorted

# Reassemble header + sorted body
cat $HEADER_FILE ${BODY_FILE}.sorted > $OUTPUT_FILE

# Clean up
rm $HEADER_FILE $BODY_FILE ${BODY_FILE}.sorted

echo "Sorted GFF3 file saved as $OUTPUT_FILE"
