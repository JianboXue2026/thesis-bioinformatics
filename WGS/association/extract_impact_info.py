"""
extract_impact_info.py
Extract HIGH / MODERATE impact variants (with gene name and variant type) from
a SnpEff-annotated VCF into two readable CSV files — used on the
PLINK-significant variant VCFs produced by extract_significant_variants.sh.

Usage:
    python extract_impact_info.py
    # edit VCF_FILE in the user-configurable section below
    # (a *_vep_annoed.py variant parses the CSQ field instead of ANN)
"""

import pysam
import csv
import os

# --- User-configurable section ----------------------------------------------
# SnpEff-annotated VCF to extract impact information from
VCF_FILE = "/path/to/snpeff_anno/all-hg38-GATK-fil-snpeff_annoed.vcf.gz"
# Output directory for the two CSV files
OUTPUT_DIR = "."
# -----------------------------------------------------------------------------


def extract_impact_info(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    high_impact_list = []
    moderate_impact_list = []
    for record in vcf:
        if not record.alts:          # no ALT (e.g. <NON_REF> records) — skip
            continue
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]
        info_dict = {key: value for key, value in record.info.items()}
        info_str = ";".join([f"{k}={v}" for k, v in info_dict.items()])
        ann_info = info_dict.get("ANN", ("N/A",))  # ANN as a tuple
        for ann_entry in ann_info:
            fields = ann_entry.split("|")
            if len(fields) < 5:
                continue
            impact = fields[2]
            gene_name = fields[3]
            variant_type = fields[1]
            if impact == "HIGH":
                high_impact_list.append(
                    [chrom, pos, ref, alt, gene_name, variant_type, impact, info_str])
            elif impact == "MODERATE":
                moderate_impact_list.append(
                    [chrom, pos, ref, alt, gene_name, variant_type, impact, info_str])
    return high_impact_list, moderate_impact_list


def save_to_csv(data, output_file):
    with open(output_file, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Chromosome", "Position", "Ref", "Alt",
                            "Gene Name", "Variant Type", "Impact", "Info"])
        csvwriter.writerows(data)
    print(f"Results have been saved to {output_file}")


def main():
    vcf_file = VCF_FILE
    base_name = os.path.basename(vcf_file)
    file_name, _ = os.path.splitext(base_name)

    high_impact_csv = os.path.join(OUTPUT_DIR, f"{file_name}_high_impact.csv")
    moderate_impact_csv = os.path.join(OUTPUT_DIR, f"{file_name}_moderate_impact.csv")

    high_impact_list, moderate_impact_list = extract_impact_info(vcf_file)
    save_to_csv(high_impact_list, high_impact_csv)
    save_to_csv(moderate_impact_list, moderate_impact_csv)


if __name__ == "__main__":
    main()
