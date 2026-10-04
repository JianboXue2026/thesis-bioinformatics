"""
extract_bamstat_info.py
Summarise `samtools stats` reports (one per sample x reference genome) into a
single CSV table for downstream mapping-quality comparison.

Expected input files (produced by samtools_index_stats.sh):
    <sample>-hg38-stat.txt / <sample>-CHM13-stat.txt / <sample>-YAO-stat.txt

Metrics are read from fixed line numbers of the samtools stats output
(samtools 1.x "SN" section layout). If your samtools version emits a
different layout, verify the line numbers before running.

Output: Data_extraction_CAP001-CAP020.csv
    Names, Ref_seq, Total_seqs, reads mapped, reads mapped and paired,
    Reads Unmapped, Reads properly paired, Percentage of Reads properly paired
"""

import csv


def extract_lines(file_path):
    """Read key metrics from one samtools stats report by fixed line numbers."""
    try:
        with open(file_path, "r", encoding="utf-8") as txt_file:
            lines = txt_file.readlines()

            total_seq_data = lines[8 - 1].strip()[24:]          # total sequences
            mapped = lines[14 - 1].strip()[17:]                 # reads mapped
            mapped_paired = lines[15 - 1].strip()[28:-52]       # reads mapped and paired
            unmapped_data = lines[16 - 1].strip()[19:]          # reads unmapped
            properly_paired_data = lines[17 - 1].strip()[26:-22]  # reads properly paired
            percentage_pro_pai_data = lines[45 - 1].strip()[44:]  # properly paired (%)

            return (total_seq_data, mapped, mapped_paired, unmapped_data,
                    properly_paired_data, percentage_pro_pai_data)
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"Error while extracting information: {e}")


def write_to_csv(file_path, column_names, data):
    try:
        with open(file_path, "w", newline="", encoding="utf-8") as csv_file:
            writer = csv.writer(csv_file, delimiter=",")
            writer.writerow(column_names)
            for row in data:
                writer.writerow(row)
            print(f"Data successfully written to {file_path}")
    except Exception as e:
        print(f"Error while writing CSV file: {e}")


# --- User-configurable section ---------------------------------------------
# Sample list: CAP001 - CAP020 (extend / shorten as needed)
name_list = ["CAP" + f"{i:03d}" for i in range(1, 20 + 1)]
# Directory holding the samtools stats reports
bamstat_dir = "../bwa-samtools-stats/"
# ---------------------------------------------------------------------------

file_list = []
data_list = []
for name in name_list:
    file_group = [bamstat_dir + name + "-CHM13-stat.txt",
                  bamstat_dir + name + "-hg38-stat.txt",
                  bamstat_dir + name + "-YAO-stat.txt"]
    file_list.append(file_group)

for file_group in file_list:
    for file in file_group:
        result = list(extract_lines(file))
        # Column 1: sample ID (first 6 chars of the file name)
        result.insert(0, file[len(bamstat_dir):len(bamstat_dir) + 6])
        # Column 2: reference genome tag (e.g. CHM13 / hg38 / YAO)
        result.insert(0, file[len(bamstat_dir) + 7:-9])
        data_list.append(result)

columns = ["Names", "Ref_seq", "Total_seqs", "reads mapped",
           "reads mapped and paired", "Reads Unmapped",
           "Reads properly paired", "Percentage of Reads properly paired"]
write_to_csv(file_path="Data_extraction_CAP001-CAP020.csv",
             column_names=columns, data=data_list)
