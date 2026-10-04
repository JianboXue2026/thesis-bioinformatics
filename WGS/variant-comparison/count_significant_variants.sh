#!/bin/bash
# count_significant_variants.sh
# Count variants with P < 0.05 in a variant-comparison output CSV
# (column 10 = P-value), e.g. produced by vcf_info_extraction_*.py.
#
# Usage: ./count_significant_variants.sh <csv_file>

csv_file=$1
count=$(awk -F',' 'NR > 1 && $10 < 0.05 {count++} END {print count}' "$csv_file")
echo "Number of entries with P-value < 0.05: ${count}"
