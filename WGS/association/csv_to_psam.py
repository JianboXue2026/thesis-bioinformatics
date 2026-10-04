"""
csv_to_psam.py
Convert a two-column sex-information CSV into a PLINK2 .psam file.

Input CSV format (with header):
    sample,sex
    CAP001,1
    CAP002,2
    ...   (sex: 1 = male, 2 = female, 0 = unknown)
"""

import csv

# --- User-configurable section ----------------------------------------------
input_csv = "/path/to/GATK-PLINK/00.trans2bed/NSCP_sex_info.csv"
output_psam = "/path/to/GATK-PLINK/00.trans2bed/NSCP_sex_info.psam"
# -----------------------------------------------------------------------------

with open(input_csv, "r") as csvfile, open(output_psam, "w") as psamfile:
    reader = csv.reader(csvfile)
    next(reader)                       # skip the CSV header
    psamfile.write("#FID IID SEX\n")   # .psam header
    for row in reader:
        sample, sex = row
        psamfile.write(f"1 {sample} {sex}\n")
