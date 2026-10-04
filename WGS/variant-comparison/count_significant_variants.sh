#!/bin/bash
# count_significant_variants.sh
# Count variants with P < 0.05 in a variant-comparison output CSV produced by
# vcf_info_extraction_*.py.
#
# The CSV is parsed with the Python csv module rather than `awk -F','` because
# the Info column may itself contain commas (e.g. multi-entry ANN/CSQ values),
# which would silently shift columns and corrupt the P-value lookup.

# Usage: ./count_significant_variants.sh <csv_file>

csv_file=$1
if [ -z "$csv_file" ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
fi

count=$(python3 - "$csv_file" <<'PY'
import csv
import sys

path = sys.argv[1]
count = 0
with open(path, newline="", encoding="utf-8") as fh:
    reader = csv.reader(fh)
    try:
        header = next(reader)
    except StopIteration:
        print(0)
        sys.exit(0)
    # Locate the P-value column by name; fall back to index 9 (0-based)
    try:
        p_idx = header.index("P-value")
    except ValueError:
        p_idx = 9
    for row in reader:
        try:
            if float(row[p_idx]) < 0.05:
                count += 1
        except (ValueError, IndexError):
            continue
print(count)
PY
)

echo "Number of entries with P-value < 0.05: ${count}"
