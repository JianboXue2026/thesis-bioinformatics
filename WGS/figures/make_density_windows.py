#!/usr/bin/env python3
"""Generate per-chromosome 1 Mb window position tables for the variant-
density ideograms (`chromosome_density_ideogram.R`).

For every primary chromosome of a reference genome the script emits one row
per 1 Mb window (Chr, Start, End) plus an empty `Count` column. Fill `Count`
with the variant counts per window (counted upstream from the filtered cohort
VCF, e.g. with `bcftools index --stats` or bedtools coverage), then normalise
to 0-1 for the `Value` column consumed by the ideogram script.

Usage
-----
    python make_density_windows.py --fasta reference.fasta \
        --output data/chromosome-density/windows.csv

The chromosome lengths are taken from a reference FASTA (.fai is created if
absent); alternatively pass --lengths-csv with columns Chr,End.
"""

import argparse
import csv
import os
import subprocess

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

WINDOW = 1_000_000  # window size in bp

# ---------------------------------------------------------------------------


def lengths_from_fasta(fasta_path):
    """Chromosome lengths from a FASTA index (building the .fai if needed)."""
    fai_path = fasta_path + ".fai"
    if not os.path.exists(fai_path):
        # List form (no shell) keeps paths with spaces safe
        subprocess.run(["samtools", "faidx", fasta_path], check=True)
    lengths = {}
    with open(fai_path, "r", encoding="utf-8") as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def lengths_from_csv(lengths_csv):
    """Chromosome lengths from a two-column CSV (Chr, End)."""
    lengths = {}
    with open(lengths_csv, "r", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lengths[row["Chr"]] = int(row["End"])
    return lengths


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--fasta", help="reference genome FASTA")
    source.add_argument("--lengths-csv", help="CSV with columns Chr,End")
    parser.add_argument("--window", type=int, default=WINDOW,
                        help="window size in bp (default: %d)" % WINDOW)
    parser.add_argument("--output", required=True, help="output CSV path")
    args = parser.parse_args(argv)

    lengths = (lengths_from_fasta(args.fasta) if args.fasta
               else lengths_from_csv(args.lengths_csv))

    with open(args.output, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Chr", "Start", "End", "Count"])
        for chrom in sorted(lengths, key=lambda c: (len(c), c)):
            start = 1
            while start <= lengths[chrom]:
                end = min(start + args.window - 1, lengths[chrom])
                writer.writerow([chrom, start, end, ""])
                start = end + 1
    print(f"Window table written: {args.output} (fill the Count column with "
          f"per-window variant counts, then normalise to 0-1 as Value)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
