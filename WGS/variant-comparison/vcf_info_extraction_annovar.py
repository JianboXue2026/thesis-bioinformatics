"""
vcf_info_extraction_annovar.py
ANNOVAR counterpart of vcf_info_extraction_snpeff.py.

ANNOVAR-annotated VCFs have no Impact field; instead they carry
ExonicFunc.refGene and AAChange.refGene. The dedicated CSVs produced here list
variants involving exonic function changes and amino-acid changes, together
with gene name and region (Func.refGene / Gene.refGene).
"""

import pysam
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- User-configurable section ------------------------------------------------
vcf_file1 = "/path/to/anno_vcf/combined_NSCP_final-YAO_annovar_anno.vcf.gz"
vcf_file2 = "/path/to/anno_vcf/combined_SCAP_final-YAO_annovar_anno.vcf.gz"
output_dir = "/path/to/results/Intra-ref_diff/"
output_csv = output_dir + "NSCP_vs_SCAP-YAO-ANNOVAR-output.csv"
unique_to_file1_vcf = output_dir + "unique_to_NSCP-YAO-ANNOVAR.vcf"
unique_to_file2_vcf = output_dir + "unique_to_SCAP-YAO-ANNOVAR.vcf"
exonic_func_csv = output_dir + "NSCP_vs_SCAP-YAO-ANNOVAR_exonic_func.csv"
aa_change_csv = output_dir + "NSCP_vs_SCAP-YAO-ANNOVAR_aa_change.csv"
required_chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
# -----------------------------------------------------------------------------


def extract_ac_an_info(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    ac_an_info = {}
    records = {}
    for record in vcf:
        if not record.alts:                     # no ALT (e.g. <NON_REF>) — skip
            continue
        chrom = record.chrom
        pos = record.pos
        ref = record.ref
        alt = record.alts[0]                    # first ALT only
        ac = (record.info.get("AC") or [0])[0]  # first AC only
        an = record.info.get("AN", 0)
        key = (chrom, pos, ref, alt)
        info_dict = {k: v for k, v in record.info.items()}
        info_str = ";".join([f"{k}={v}" for k, v in info_dict.items()])
        ann_info = info_dict.get("ANN", ("N/A",))
        exonic_func = info_dict.get("ExonicFunc.refGene", ".")[0]
        aa_change = info_dict.get("AAChange.refGene", ".")[0]
        ac_an_info[key] = {
            "AC": ac, "AN": an, "INFO": info_str, "ANN": ann_info,
            "FUNC": info_dict.get("Func.refGene", ".")[0],
            "GENE": info_dict.get("Gene.refGene", ".")[0],
            "EXONICFUNC": exonic_func,
            "AACHANGE": aa_change,
        }
        records[key] = record
    return ac_an_info, records


ac_an_info1, records1 = extract_ac_an_info(vcf_file1)
ac_an_info2, records2 = extract_ac_an_info(vcf_file2)

all_keys = set(ac_an_info1.keys()) | set(ac_an_info2.keys())


def contingency_test(table):
    """Chi-square test, falling back to Fisher's exact test when the smallest
    EXPECTED cell count is < 5. Returns (chi2, p_value, test_type)."""
    total = int(table.sum())
    if total == 0:
        return "N/A", float("nan"), "Skipped (empty table)"
    expected = (table.sum(axis=1, keepdims=True)
                * table.sum(axis=0, keepdims=True) / total)
    if expected.min() < 5:
        _, p = fisher_exact(table)
        return "N/A", p, "Fisher Exact"
    chi2, p, _, _ = chi2_contingency(table)
    return chi2, p, "Chi-Square"


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
    func = (ac_an_info1[key]["FUNC"] if key in ac_an_info1
            else ac_an_info2[key]["FUNC"])
    gene = (ac_an_info1[key]["GENE"] if key in ac_an_info1
            else ac_an_info2[key]["GENE"])
    exonic_func = (ac_an_info1[key]["EXONICFUNC"] if key in ac_an_info1
                   else ac_an_info2[key]["EXONICFUNC"])
    aa_change = (ac_an_info1[key]["AACHANGE"] if key in ac_an_info1
                 else ac_an_info2[key]["AACHANGE"])

    # Strip the surrounding tuple quotes pysam returns for these fields
    func = func.strip("('").strip("',)")
    gene = gene.strip("('").strip("',)")
    exonic_func = exonic_func.strip("('").strip("',)")
    aa_change = aa_change.strip("('").strip("',)")

    table = np.array([[ac1, an1 - ac1], [ac2, an2 - ac2]])
    chi2, p, test_type = contingency_test(table)

    combined_result = [chrom, pos, ref, alt, ac1, an1, ac2, an2, chi2, p,
                       ";".join(ann), info, test_type]

    exonic_func_results = []
    aa_change_results = []
    if exonic_func and exonic_func != ".":
        exonic_func_results.append(
            (chrom, pos, ref, alt, ac1, an1, ac2, an2, chi2, p, info,
             test_type, gene, func, exonic_func))
    if aa_change and aa_change != ".":
        aa_change_results.append(
            (chrom, pos, ref, alt, ac1, an1, ac2, an2, chi2, p, info,
             test_type, gene, func, aa_change))

    return combined_result, exonic_func_results, aa_change_results, chrom


results = []
exonic_func_results = []
aa_change_results = []
chromosome_data = defaultdict(list)

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process_variant, key): key for key in all_keys}
    for future in as_completed(futures):
        combined_result, exonic_func, aa_change, chrom = future.result()
        results.append(combined_result)
        exonic_func_results.extend(exonic_func)
        aa_change_results.extend(aa_change)
        chromosome_data[chrom].append(combined_result)

results.sort(key=lambda x: (x[0], x[1]))
exonic_func_results.sort(key=lambda x: (x[0], x[1]))
aa_change_results.sort(key=lambda x: (x[0], x[1]))

header_all = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
              "Chi2", "P-value", "ANN", "Info", "Test Type"]
header_detail = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
                 "Chi2", "P-value", "Info", "Test Type", "Gene Name", "Variant Type"]

with open(output_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_all)
    csvwriter.writerows(results)
print(f"Results have been saved to {output_csv}")

with open(exonic_func_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_detail + ["ExonicFunc"])
    csvwriter.writerows(exonic_func_results)
print(f"Exonic function results have been saved to {exonic_func_csv}")

with open(aa_change_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_detail + ["AAChange"])
    csvwriter.writerows(aa_change_results)
print(f"AA change results have been saved to {aa_change_csv}")

for chrom in required_chromosomes:
    if chrom in chromosome_data:
        data = chromosome_data[chrom]
        data.sort(key=lambda x: x[1])
        chrom_csv = output_dir + f"NSCP_vs_SCAP-YAO-ANNOVAR-output_{chrom}.csv"
        with open(chrom_csv, "w", newline="") as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(header_all)
            csvwriter.writerows(data)
        print(f"Chromosome {chrom} results have been saved to {chrom_csv}")

# Unique variants -> VCF
unique_to_file1 = set(ac_an_info1.keys()) - set(ac_an_info2.keys())
vcf1 = pysam.VariantFile(vcf_file1)
with pysam.VariantFile(unique_to_file1_vcf, "w", header=vcf1.header) as out_vcf1:
    for key in unique_to_file1:
        out_vcf1.write(records1[key])
print(f"Unique variants in {vcf_file1} have been saved to {unique_to_file1_vcf}")

unique_to_file2 = set(ac_an_info2.keys()) - set(ac_an_info1.keys())
vcf2 = pysam.VariantFile(vcf_file2)
with pysam.VariantFile(unique_to_file2_vcf, "w", header=vcf2.header) as out_vcf2:
    for key in unique_to_file2:
        out_vcf2.write(records2[key])
print(f"Unique variants in {vcf_file2} have been saved to {unique_to_file2_vcf}")
