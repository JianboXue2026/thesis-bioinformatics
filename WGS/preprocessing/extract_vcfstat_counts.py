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

            lines_remainder = txt_lines_count // interval
            line_num_list_snp = [start_value_snp + i * interval for i in range(lines_remainder)]
            line_num_list_indel = [start_value_indel + i * interval for i in range(lines_remainder)]

            snp_list = []
            indel_list = []

            for line_num_snp in line_num_list_snp:
                snp_list.append(lines[line_num_snp - 1].strip()[22:29])
            for line_num_indel in line_num_list_indel:
                indel_list.append(lines[line_num_indel - 1].strip()[22:28])

            return snp_list, indel_list

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"Error while extracting information: {e}")


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
