"""
vcf_info_extraction_snpeff.py
Compare two SnpEff-annotated cohort VCFs (e.g. NSCP vs SCAP, same reference
genome) variant by variant:

  1. Extract AC / AN / INFO / ANN for every variant in both VCFs.
  2. For each variant build a 2x2 contingency table from allele counts:
         [ AC_group1, AN_group1 - AC_group1 ]
         [ AC_group2, AN_group2 - AC_group2 ]
     Variants present in only one group are kept as well (missing group gets
     AC = 0, AN = the other group's AN).
  3. Test each table with a chi-square test; switch to Fisher's exact test
     when any observed count in the 2x2 table is < 5.
  4. Write:
       - one combined CSV (all variants, sorted by chromosome & position)
       - per-chromosome CSVs (chr1-chr22, chrX, chrY, chrM only)
       - HIGH / MODERATE impact CSVs with gene name and variant type parsed
         from the SnpEff ANN field
       - two VCFs containing variants unique to either group

Variants are processed in parallel with a ThreadPoolExecutor.
"""

import pysam
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- User-configurable section ------------------------------------------------
vcf_file1 = "/path/to/anno_vcf/combined_NSCP_final-YAO_anno_clean.vcf.gz"
vcf_file2 = "/path/to/anno_vcf/combined_SCAP_final-YAO_anno_clean.vcf.gz"
output_dir = "/path/to/results/Intra-ref_diff/"
output_csv = output_dir + "NSCP_vs_SCAP-YAO-output.csv"
unique_to_file1_vcf = output_dir + "unique_to_NSCP-YAO.vcf"
unique_to_file2_vcf = output_dir + "unique_to_SCAP-YAO.vcf"
high_impact_csv = output_dir + "NSCP_vs_SCAP-YAO_high_impact.csv"
moderate_impact_csv = output_dir + "NSCP_vs_SCAP-YAO_moderate_impact.csv"
# Chromosomes for the per-chromosome CSVs (restrict the huge hg38 contig list)
required_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
# -----------------------------------------------------------------------------


def extract_ac_an_info(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    ac_an_info = {}
    records = {}
    for record in vcf:
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]
        ac = record.info.get("AC", [0])[0]
        an = record.info.get("AN", 0)
        key = (chrom, pos, ref, alt)
        info_dict = {k: v for k, v in record.info.items()}
        info_str = ";".join([f"{k}={v}" for k, v in info_dict.items()])
        ann_info = info_dict.get("ANN", ("N/A",))  # ANN as a tuple
        ac_an_info[key] = {"AC": ac, "AN": an, "INFO": info_str, "ANN": ann_info}
        records[key] = record
    return ac_an_info, records


ac_an_info1, records1 = extract_ac_an_info(vcf_file1)
ac_an_info2, records2 = extract_ac_an_info(vcf_file2)

# Union of all variants from both groups (group-specific variants get AC=0)
all_keys = set(ac_an_info1.keys()) | set(ac_an_info2.keys())


def process_variant(key):
    chrom, pos, ref, alt = key
    ac1 = ac_an_info1[key]["AC"] if key in ac_an_info1 else 0
    an1 = ac_an_info1[key]["AN"] if key in ac_an_info1 else ac_an_info2[key]["AN"]
    ac2 = ac_an_info2[key]["AC"] if key in ac_an_info2 else 0
    an2 = ac_an_info2[key]["AN"] if key in ac_an_info2 else ac_an_info1[key]["AN"]
    info = (ac_an_info1[key]["INFO"] if key in ac_an_info1
            else ac_an_info2[key]["INFO"])
    ann = (ac_an_info1[key]["ANN"] if key in ac_an_info1
           else ac_an_info2[key]["ANN"])

    # 2x2 contingency table
    table = np.array([[ac1, an1 - ac1], [ac2, an2 - ac2]])

    # Chi-square, or Fisher's exact test when any observed count in the table < 5
    if np.any(table < 5):
        oddsratio, p = fisher_exact(table)
        test_type = "Fisher Exact"
        chi2 = "N/A"
    else:
        chi2, p, _, expected = chi2_contingency(table)
        test_type = "Chi-Square"

    combined_result = [chrom, pos, ref, alt, ac1, an1, ac2, an2,
                       chi2, p, ";".join(ann), info, test_type]

    # Classify by SnpEff-predicted impact; parse gene name & variant type
    high_impact = []
    moderate_impact = []
    for ann_entry in ann:
        fields = ann_entry.split("|")
        if len(fields) < 5:
            continue
        impact = fields[2]
        gene_name = fields[3]
        variant_type = fields[1]
        detailed_result = (chrom, pos, ref, alt, ac1, an1, ac2, an2,
                           chi2, p, info, test_type, gene_name, variant_type, impact)
        if impact == "HIGH":
            high_impact.append(detailed_result)
        elif impact == "MODERATE":
            moderate_impact.append(detailed_result)

    return combined_result, high_impact, moderate_impact, chrom


results = []
high_impact_results = []
moderate_impact_results = []
chromosome_data = defaultdict(list)

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process_variant, key): key for key in all_keys}
    for future in as_completed(futures):
        combined_result, high_impact, moderate_impact, chrom = future.result()
        results.append(combined_result)
        high_impact_results.extend(high_impact)
        moderate_impact_results.extend(moderate_impact)
        chromosome_data[chrom].append(combined_result)

# Sort by chromosome, then position
results.sort(key=lambda x: (x[0], x[1]))
high_impact_results.sort(key=lambda x: (x[0], x[1]))
moderate_impact_results.sort(key=lambda x: (x[0], x[1]))

header_all = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
              "Chi2", "P-value", "ANN", "Info", "Test Type"]
header_impact = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
                 "Chi2", "P-value", "Info", "Test Type", "Gene Name", "Variant Type", "Impact"]

with open(output_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_all)
    csvwriter.writerows(results)
print(f"Results have been saved to {output_csv}")

with open(high_impact_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_impact)
    csvwriter.writerows(high_impact_results)
print(f"High impact results have been saved to {high_impact_csv}")

with open(moderate_impact_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_impact)
    csvwriter.writerows(moderate_impact_results)
print(f"Moderate impact results have been saved to {moderate_impact_csv}")

# Per-chromosome CSVs, restricted to the primary chromosomes
for chrom in required_chromosomes:
    if chrom in chromosome_data:
        data = chromosome_data[chrom]
        data.sort(key=lambda x: x[1])
        chrom_csv = output_dir + f"NSCP_vs_SCAP-YAO-output_{chrom}.csv"
        with open(chrom_csv, "w", newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(header_all)
            csvwriter.writerows(data)
        print(f"Chromosome {chrom} results have been saved to {chrom_csv}")

# Variants unique to group 1 -> VCF
unique_to_file1 = set(ac_an_info1.keys()) - set(ac_an_info2.keys())
vcf1 = pysam.VariantFile(vcf_file1)
with pysam.VariantFile(unique_to_file1_vcf, "w", header=vcf1.header) as out_vcf1:
    for key in unique_to_file1:
        out_vcf1.write(records1[key])
print(f"Unique variants in {vcf_file1} have been saved to {unique_to_file1_vcf}")

# Variants unique to group 2 -> VCF
unique_to_file2 = set(ac_an_info2.keys()) - set(ac_an_info1.keys())
vcf2 = pysam.VariantFile(vcf_file2)
with pysam.VariantFile(unique_to_file2_vcf, "w", header=vcf2.header) as out_vcf2:
    for key in unique_to_file2:
        out_vcf2.write(records2[key])
print(f"Unique variants in {vcf_file2} have been saved to {unique_to_file2_vcf}")
