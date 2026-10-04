"""
extract_chromosome_lengths.py
Extract chromosome lengths from a reference genome FASTA into a CSV table —
input for chromosome-level visualisations (karyoploteR / RIdeogram).

Run inside the conda environment that provides Biopython.
"""

import csv
from Bio import SeqIO


def extract_chromosome_lengths(fasta_file, output_csv):
    chromosome_lengths = {}
    # Parse lengths from the FASTA records
    for record in SeqIO.parse(fasta_file, "fasta"):
        chromosome_lengths[record.id] = len(record.seq)

    # Write to CSV
    with open(output_csv, "w", newline="") as csvfile:
        fieldnames = ["Chromosome", "Length"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for chromosome, length in chromosome_lengths.items():
            writer.writerow({"Chromosome": chromosome, "Length": length})


# --- User-configurable section ----------------------------------------------
fasta_file = "/path/to/ref_seq/YAO/T2T-Yao-hp.v1.1.fasta"
output_csv = "/path/to/results/chr_length/YAO.length.csv"
# -----------------------------------------------------------------------------

extract_chromosome_lengths(fasta_file, output_csv)
