#!/usr/bin/env python3
"""Box-line plots: alignment / QC metrics across the three reference genomes.

For each metric the script reads a wide CSV (Sample_ID, CHM13, GRCh38, YAO),
draws a box plot per reference genome, overlays one grey line per sample to
visualise the paired nature of the comparison (every sample was aligned
against all three references), annotates Q25 / median / Q75, and marks
pairwise Wilcoxon signed-rank tests with significance stars.

Input CSVs are expected in ``data/boxline/`` (example datasets shipped with
this repository; substitute your own files with identical column layout).

Usage
-----
    python box_line_plots.py                    # all metrics
    python box_line_plots.py --metric error_rate bases_mapped_cigar

Edit the user-configurable section below before running.
"""

import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from itertools import combinations
from scipy.stats import wilcoxon

# ---------------------------------------------------------------------------
# User-configurable section
# ---------------------------------------------------------------------------

DATA_DIR = "data/boxline"       # input CSVs
OUTPUT_DIR = "output/boxline"   # figures are written here

GROUPS = ["CHM13", "GRCh38", "YAO"]
PALETTE = {"CHM13": "#1f77b4", "GRCh38": "#2ca02c", "YAO": "#ff7f0e"}

# Per-metric configuration:
#   file   input CSV name inside DATA_DIR
#   ylabel y-axis label of the figure
#   qfmt   "sci" -> scientific notation (Q25/Med/Q75), "pct" -> percentage
#   ymax   anchor for the annotation block: "q95" / "q97" (quantile) or "max"
METRICS = {
    "bases_mapped_cigar": dict(file="box-line_bases_mapped_cigar.csv",
                               ylabel="Bases Mapped (bp)", qfmt="sci", ymax="q95"),
    "mismatches":         dict(file="box-line_mismatches.csv",
                               ylabel="Mismatched Bases (bp)", qfmt="sci", ymax="q97"),
    "error_rate":         dict(file="box-line_error_rate.csv",
                               ylabel="Error Rate (%)", qfmt="pct", ymax="max"),
    "mapped_and_paired":  dict(file="box-line_mapped_and_paired.csv",
                               ylabel="Properly Paired Reads", qfmt="sci", ymax="q95"),
    "reads_unmapped":     dict(file="box-line_reads_unmapped.csv",
                               ylabel="Unmapped Reads", qfmt="sci", ymax="q95"),
    "Het_Mutant":         dict(file="box-line_Het_Mutant.csv",
                               ylabel="Heterozygous Variants", qfmt="sci", ymax="q95"),
    "Hom_Mutant":         dict(file="box-line_Hom_Mutant.csv",
                               ylabel="Homozygous Variants", qfmt="sci", ymax="q95"),
    "All_Mutant":         dict(file="box-line_All_Mutant.csv",
                               ylabel="All Variants", qfmt="sci", ymax="q95"),
}

# ---------------------------------------------------------------------------


def format_p_value(p):
    """Map a p-value onto */**/*** significance stars."""
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    return 'n.s.'


def format_quantile(value, qfmt):
    """Q25/Med/Q75 annotation format: scientific or percentage."""
    if qfmt == "pct":
        return f"{value:.3f}%"
    return f"{value:.4e}"


def ymax_anchor(series, strategy):
    if strategy == "max":
        return series.max()
    quantile = {"q95": 0.95, "q97": 0.97}[strategy]
    return series.quantile(quantile)


def plot_metric(metric, cfg, data_dir=DATA_DIR, output_dir=OUTPUT_DIR):
    file_path = os.path.join(data_dir, cfg["file"])
    data = pd.read_csv(file_path)
    data_filtered = data[["Sample_ID"] + GROUPS]

    # Long format for the box plot
    data_long = data_filtered.melt(id_vars="Sample_ID", var_name="Group", value_name="Value")

    # Keep only samples whose values for every reference fall inside the
    # group-wise 5%-95% range; outliers of a single reference drop the whole
    # sample from the line overlay (the box plot itself uses the full data).
    q5 = data_long.groupby("Group")["Value"].quantile(0.05)
    q95 = data_long.groupby("Group")["Value"].quantile(0.95)
    data_long["Keep"] = (data_long["Value"] >= data_long["Group"].map(q5)) & \
                        (data_long["Value"] <= data_long["Group"].map(q95))
    valid_samples = data_long.groupby("Sample_ID")["Keep"].all()
    filtered = data_long[data_long["Sample_ID"].isin(valid_samples[valid_samples].index)] \
        .drop(columns="Keep")

    # Pairwise Wilcoxon signed-rank tests on the full (unfiltered) data
    p_values = {}
    for g1, g2 in combinations(GROUPS, 2):
        stat, p = wilcoxon(data_filtered[g1], data_filtered[g2])
        p_values[(g1, g2)] = p

    median_values = data_filtered[GROUPS].median()
    q1_values = data_filtered[GROUPS].quantile(0.25)
    q3_values = data_filtered[GROUPS].quantile(0.75)

    plt.figure(figsize=(7, 7))
    # hue="Group" + legend=False avoids the deprecation warning raised when
    # passing `palette` without `hue` (matplotlib >= 3.6 / seaborn >= 0.13)
    sns.boxplot(data=data_long, x="Group", y="Value", hue="Group",
                palette=PALETTE, legend=False, width=0.5, showfliers=False)

    y_max = ymax_anchor(data_long["Value"], cfg["ymax"])
    y_offset = (y_max - data_long["Value"].min()) * 0.3
    y_label_offset = y_max + y_offset * 0.5

    for i, group in enumerate(GROUPS):
        plt.text(i, y_label_offset - y_offset * 0.3,
                 f"Q25: {format_quantile(q1_values[group], cfg['qfmt'])}",
                 ha='center', va='bottom', color='darkred', fontsize=8, fontweight='bold')
        plt.text(i, y_label_offset,
                 f"Med: {format_quantile(median_values[group], cfg['qfmt'])}",
                 ha='center', va='bottom', color='darkred', fontsize=8, fontweight='bold')
        plt.text(i, y_label_offset + y_offset * 0.3,
                 f"Q75: {format_quantile(q3_values[group], cfg['qfmt'])}",
                 ha='center', va='bottom', color='darkred', fontsize=8, fontweight='bold')

    # One grey line per sample across the three references (5-95% kept set)
    for sample_id in filtered["Sample_ID"].unique():
        sample_data = filtered[filtered["Sample_ID"] == sample_id]
        plt.plot(sample_data["Group"], sample_data["Value"],
                 marker='o', color="black", alpha=0.5, linestyle="-")

    # Significance brackets
    for i, ((g1, g2), p) in enumerate(p_values.items()):
        y_position = y_max + y_offset + (0.5 * i) * y_offset
        plt.plot([GROUPS.index(g1), GROUPS.index(g2)], [y_position, y_position],
                 color="black")
        plt.text((GROUPS.index(g1) + GROUPS.index(g2)) / 2, y_position,
                 format_p_value(p), ha='center', va='bottom')

    plt.legend([], [], frameon=False)
    plt.xlabel("Reference Genome")
    plt.ylabel(cfg["ylabel"])

    os.makedirs(output_dir, exist_ok=True)
    save_path = os.path.join(output_dir, f"box-line_{metric}.pdf")
    plt.savefig(save_path, dpi=300, bbox_inches='tight', format='pdf')
    plt.close()
    print(f"Saved: {save_path}")
    return p_values


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
        plot_metric(metric, METRICS[metric])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
