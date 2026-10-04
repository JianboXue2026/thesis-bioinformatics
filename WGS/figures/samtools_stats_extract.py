#!/usr/bin/env python3
"""Extract per-sample metrics from samtools stats reports.

For every ``<sample_id>_<reference>.stats`` file in a directory, pull one
selected ``SN`` line (e.g. "reads mapped and paired") and collect the values
into a long-format CSV (sample, reference, value). An optional wide-format
CSV (one column per reference) can be written for the box-line plots in
``box_line_plots.py``.

Usage
-----
    python samtools_stats_extract.py \
        --stats-dir  /path/to/01.Samtools_stats \
        --metric     mapped_and_paired \
        --output     data/boxline/raw_mapped_and_paired.csv \
        --wide-output data/boxline/box-line_mapped_and_paired.csv

Edit the user-configurable section below before running.
"""

import argparse
import csv
import os
import re
import sys

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

# Directory that holds one samtools stats report per sample x reference,
# named "<sample_id>_<reference>.stats" (e.g. "CAP001_CHM13.stats").
DEFAULT_STATS_DIR = "/path/to/01.Samtools_stats"

# Default output locations (relative to the figures/ directory).
DEFAULT_OUTPUT = "data/boxline/raw_<metric>.csv"
DEFAULT_WIDE_OUTPUT = "data/boxline/box-line_<metric>.csv"

# SN metric keyword (must match the samtools stats wording exactly) -> short key
SN_METRICS = {
    "reads mapped and paired": "mapped_and_paired",
    "reads unmapped": "reads_unmapped",
    "bases mapped (cigar)": "bases_mapped_cigar",
    "mismatches": "mismatches",
    "error rate": "error_rate",
    "raw total sequences": "total_reads",
    "reads mapped": "reads_mapped",
    "average quality": "average_quality",
}

# Reference label order for the wide output. References not listed here are
# still exported in the long CSV but excluded from the wide table.
REF_ORDER = ["CHM13", "GRCh38", "YAO"]

# ---------------------------------------------------------------------------


def parse_file_name(file_name):
    """Return (sample_id, reference) from '<sample>_<ref>.stats'.

    Raises ValueError when the name does not follow the convention so the
    problem surfaces immediately instead of producing silently wrong labels.
    """
    m = re.match(r"^(?P<sample>.+)_(?P<ref>.+)\.stats$", file_name)
    if not m:
        raise ValueError(
            f"File name '{file_name}' does not match '<sample>_<ref>.stats'"
        )
    return m.group("sample"), m.group("ref")


def extract_sn_value(stats_path, keyword):
    """Return the numeric value of the SN line containing ``keyword``.

    Values are parsed as float; callers can re-format as needed. Returns
    None when the line is absent (e.g. empty BAM).
    """
    with open(stats_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("SN") and keyword in line:
                return float(line.rstrip("\n").split("\t")[2])
    return None


def collect(stats_dir, keyword):
    """Walk the stats directory and collect (sample, ref, value) triples."""
    rows = []
    for file_name in sorted(os.listdir(stats_dir)):
        if not file_name.endswith(".stats"):
            continue
        sample_id, reference = parse_file_name(file_name)
        value = extract_sn_value(os.path.join(stats_dir, file_name), keyword)
        rows.append((sample_id, reference, value))
    return rows


def save_long(rows, output_csv):
    with open(output_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Sample_ID", "Reference", "Value"])
        writer.writerows(rows)


def save_wide(rows, output_csv):
    """Pivot the long rows into one column per reference."""
    columns = list(REF_ORDER)
    by_sample = {}
    for sample_id, reference, value in rows:
        if reference not in columns:
            continue
        by_sample.setdefault(sample_id, {})[reference] = value
    with open(output_csv, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Sample_ID"] + columns)
        for sample_id in sorted(by_sample):
            writer.writerow(
                [sample_id] + [by_sample[sample_id].get(r, "") for r in columns]
            )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--stats-dir", default=DEFAULT_STATS_DIR,
                        help="directory with <sample>_<ref>.stats files")
    parser.add_argument("--metric", choices=sorted(set(SN_METRICS.values())),
                        help="which SN metric to extract (required unless "
                             "--list-metrics is given)")
    parser.add_argument("--output", default=None,
                        help="long-format CSV output path (default: %s)" % DEFAULT_OUTPUT)
    parser.add_argument("--wide-output", default=None,
                        help="optional wide-format CSV output path")
    parser.add_argument("--list-metrics", action="store_true",
                        help="print available metrics and exit")
    args = parser.parse_args(argv)

    if args.list_metrics:
        for keyword, key in sorted(SN_METRICS.items(), key=lambda kv: kv[1]):
            print(f"{key:22s} -> SN line containing '{keyword}'")
        return 0

    if args.metric is None:
        parser.error("--metric is required unless --list-metrics is given")

    keyword = next(k for k, v in SN_METRICS.items() if v == args.metric)
    output = args.output or DEFAULT_OUTPUT.replace("<metric>", args.metric)

    rows = collect(args.stats_dir, keyword)
    if not rows:
        sys.exit(f"No .stats files found in {args.stats_dir}")
    save_long(rows, output)
    print(f"Long CSV written:   {output}  ({len(rows)} rows)")

    if args.wide_output is None:
        args.wide_output = DEFAULT_WIDE_OUTPUT.replace("<metric>", args.metric)
    save_wide(rows, args.wide_output)
    print(f"Wide CSV written:   {args.wide_output}")

    missing = sum(1 for _, _, v in rows if v is None)
    if missing:
        print(f"WARNING: {missing} stats files did not contain '{keyword}'.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
