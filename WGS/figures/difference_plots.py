#!/usr/bin/env python3
"""Per-sample difference scatter plots: T2T-YAO vs the other two references.

For each metric the script reads a CSV with columns

    Sample_ID, YAO-CHM13 (Mb), YAO-GRCh38 (Mb)

i.e. how much smaller/larger the metric is on T2T-YAO than on T2T-CHM13 or
GRCh38 (YAO as the baseline). Samples are sorted by the YAO-GRCh38 column so
the two series can be compared sample by sample, and a scatter of the two
differences is drawn per sample.

Input CSVs are expected in ``data/difference/`` (example datasets shipped
with this repository; substitute your own files with identical layout).

Usage
-----
    python difference_plots.py                      # all three metrics
    python difference_plots.py --metric mapped_bases

Edit the user-configurable section below before running.
"""

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

DATA_DIR = "data/difference"
OUTPUT_DIR = "output/difference"

# key -> (input CSV, y-axis label of the PNG, y-axis label of the renewed PDF)
METRICS = {
    "mapped_bases": dict(
        file="difference_mapped_bases.csv",
        ylabel_png="Base Reduction of Aligned Reads (Mb)",
        ylabel_pdf="Bases Reduction of Aligned (Mb)"),
    "mismatches": dict(
        file="difference_mismatches.csv",
        ylabel_png="Base Mismatch Increase of Aligned Reads (Mb)",
        ylabel_pdf="Bases Mismatch Increase of Aligned (Mb)"),
    "percentile_error_rate": dict(
        file="difference-percentile_error_rate.csv",
        ylabel_png="Increase Percentile of Error Rate (%)",
        ylabel_pdf="Increase Percentile of Error Rate (%)"),
}

# ---------------------------------------------------------------------------


def plot_difference(metric, cfg, data_dir=DATA_DIR, output_dir=OUTPUT_DIR):
    file_path = os.path.join(data_dir, cfg["file"])
    df = pd.read_csv(file_path)

    df_sorted = df.sort_values(by="YAO-GRCh38 (Mb)", ascending=True)
    sample_ids = df_sorted["Sample_ID"]
    chm13_values = df_sorted["YAO-CHM13 (Mb)"]
    grch38_values = df_sorted["YAO-GRCh38 (Mb)"]

    outputs = []
    for suffix, ylabel, fmt in (("png", cfg["ylabel_png"], "png"),
                                ("renew", cfg["ylabel_pdf"], "pdf")):
        plt.figure(figsize=(18, 6))
        if suffix == "png":
            # Original style: show both differences plus the YAO baseline at 0
            plt.scatter(sample_ids, chm13_values, label="CHM13", color="#B8860B", s=10)
            plt.scatter(sample_ids, grch38_values, label="GRCh38", color="#00008B", s=10)
            plt.scatter(sample_ids, [0] * len(sample_ids),
                        color="#8B0000", label="YAO (Baseline)", s=10)
            plt.gca().invert_yaxis()
        else:
            # Renewed style: reduction of CHM13 / GRCh38 relative to YAO only
            plt.scatter(sample_ids, chm13_values, label="CHM13-YAO", color="#B8860B", s=10)
            plt.scatter(sample_ids, grch38_values, label="GRCh38-YAO", color="#00008B", s=10)

        plt.xticks(rotation=90, fontsize=8)
        plt.xlabel("Sample_ID", fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        plt.legend(fontsize=14)
        plt.tight_layout()

        os.makedirs(output_dir, exist_ok=True)
        save_path = os.path.join(output_dir, f"difference_{metric}-{suffix}.{fmt}")
        plt.savefig(save_path, dpi=300, format=fmt)
        plt.close()
        outputs.append(save_path)

    for path in outputs:
        print(f"Saved: {path}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--metric", nargs="*", default=None,
                        help="metric keys to plot (default: all); "
                             "available: " + ", ".join(METRICS))
    args = parser.parse_args(argv)

    selected = args.metric if args.metric else list(METRICS)
    unknown = [m for m in selected if m not in METRICS]
    if unknown:
        parser.error("unknown metric(s): " + ", ".join(unknown))

    for metric in selected:
        plot_difference(metric, METRICS[metric])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
