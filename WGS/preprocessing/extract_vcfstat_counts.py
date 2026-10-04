"""
extract_vcfstat_counts.py
Extract SNP and INDEL counts from a bcftools-stats-style log file
(generated from per-chromosome / per-batch VCF statistics) and write them
into a two-column CSV.

The counts are read from fixed line numbers with a fixed interval between
consecutive records — verify both against your own stats file before running.

Input : filtered_stat-01.log.txt (edit the path below)
Output: snp_indel_counts_output.csv  (SNP count, INDEL count)
"""

import csv


def extract_lines(file_path):
    try:
        with open(file_path, "r", encoding="utf-8") as txt_file:
            lines = txt_file.readlines()

        txt_lines_count = len(lines)
        start_value_snp = 13    # first line holding a SNP count
        start_value_indel = 26  # first line holding an INDEL count
        interval = 26           # line interval between consecutive records

        # Enumerate only line numbers that actually exist, instead of using an
        # integer division that can overshoot the end of the file.
        line_num_list_snp = list(range(start_value_snp, txt_lines_count + 1, interval))
        line_num_list_indel = list(range(start_value_indel, txt_lines_count + 1, interval))

        snp_list = [lines[n - 1].strip()[22:29] for n in line_num_list_snp]
        indel_list = [lines[n - 1].strip()[22:28] for n in line_num_list_indel]

        if not snp_list:
            print(f"Warning: no SNP records found in '{file_path}' — check the layout.")
        return snp_list, indel_list

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
        return [], []
    except Exception as e:
        print(f"Error while extracting information: {e}")
        return [], []


# --- User-configurable section ---------------------------------------------
snp_list, indel_list = extract_lines(file_path="filtered_stat-01.log.txt")
# ---------------------------------------------------------------------------

csv_filename = "snp_indel_counts_output.csv"

with open(csv_filename, "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(["SNP_count", "INDEL_count"])
    for item1, item2 in zip(snp_list, indel_list):
        csv_writer.writerow([item1, item2])

print(f"Data has been written to {csv_filename}.")
