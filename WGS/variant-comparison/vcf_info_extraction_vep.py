"""
vcf_info_extraction_vep.py
VEP counterpart of vcf_info_extraction_snpeff.py.

VEP-annotated VCFs store annotations in the CSQ field, whose layout is similar
to SnpEff's ANN (pipe-separated). VEP additionally reports a MODIFIER impact
tier, which is written to its own CSV alongside HIGH / MODERATE / LOW.
"""

import pysam
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
import csv
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- User-configurable section ------------------------------------------------
vcf_file1 = "/path/to/anno_vcf/combined_NSCP_final-YAO_vep_anno.vcf.gz"
vcf_file2 = "/path/to/anno_vcf/combined_SCAP_final-YAO_vep_anno.vcf.gz"
output_dir = "/path/to/results/Intra-ref_diff/"
output_csv = output_dir + "NSCP_vs_SCAP-YAO-vep-output.csv"
unique_to_file1_vcf = output_dir + "unique_to_NSCP-YAO-vep.vcf"
unique_to_file2_vcf = output_dir + "unique_to_SCAP-YAO-vep.vcf"
high_impact_csv = output_dir + "NSCP_vs_SCAP-YAO-vep_high_impact.csv"
moderate_impact_csv = output_dir + "NSCP_vs_SCAP-YAO-vep_moderate_impact.csv"
low_impact_csv = output_dir + "NSCP_vs_SCAP-YAO-vep_low_impact.csv"
modifier_impact_csv = output_dir + "NSCP_vs_SCAP-YAO-vep_modifier_impact.csv"
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
        csq_info = info_dict.get("CSQ", ("N/A",))  # CSQ as a tuple
        ac_an_info[key] = {"AC": ac, "AN": an, "INFO": info_str, "CSQ": csq_info}
        records[key] = record
    return ac_an_info, records


ac_an_info1, records1 = extract_ac_an_info(vcf_file1)
ac_an_info2, records2 = extract_ac_an_info(vcf_file2)

all_keys = set(ac_an_info1.keys()) | set(ac_an_info2.keys())


def process_variant(key):
    chrom, pos, ref, alt = key
    ac1 = ac_an_info1[key]["AC"] if key in ac_an_info1 else 0
    an1 = ac_an_info1[key]["AN"] if key in ac_an_info1 else ac_an_info2[key]["AN"]
    ac2 = ac_an_info2[key]["AC"] if key in ac_an_info2 else 0
    an2 = ac_an_info2[key]["AN"] if key in ac_an_info2 else ac_an_info1[key]["AN"]
    info = (ac_an_info1[key]["INFO"] if key in ac_an_info1
            else ac_an_info2[key]["INFO"])
    csq = (ac_an_info1[key]["CSQ"] if key in ac_an_info1
           else ac_an_info2[key]["CSQ"])

    table = np.array([[ac1, an1 - ac1], [ac2, an2 - ac2]])
    if np.any(table < 5):
        oddsratio, p = fisher_exact(table)
        test_type = "Fisher Exact"
        chi2 = "N/A"
    else:
        chi2, p, _, expected = chi2_contingency(table)
        test_type = "Chi-Square"

    combined_result = [chrom, pos, ref, alt, ac1, an1, ac2, an2,
                       chi2, p, ";".join(csq), info, test_type]

    high_impact, moderate_impact, low_impact, modifier_impact = [], [], [], []
    for csq_entry in csq:
        fields = csq_entry.split("|")
        if len(fields) < 5:
            continue
        impact = fields[2]
        gene_name = fields[3]
        variant_type = fields[1]
        detailed_result = (chrom, pos, ref, alt, ac1, an1, ac2, an2, chi2, p,
                           info, test_type, gene_name, variant_type, impact)
        if impact == "HIGH":
            high_impact.append(detailed_result)
        elif impact == "MODERATE":
            moderate_impact.append(detailed_result)
        elif impact == "LOW":
            low_impact.append(detailed_result)
        elif impact == "MODIFIER":
            modifier_impact.append(detailed_result)

    return (combined_result, high_impact, moderate_impact, low_impact,
            modifier_impact, chrom)


results = []
high_impact_results = []
moderate_impact_results = []
low_impact_results = []
modifier_impact_results = []
chromosome_data = defaultdict(list)

with ThreadPoolExecutor() as executor:
    futures = {executor.submit(process_variant, key): key for key in all_keys}
    for future in as_completed(futures):
        (combined_result, high_impact, moderate_impact, low_impact,
         modifier_impact, chrom) = future.result()
        results.append(combined_result)
        high_impact_results.extend(high_impact)
        moderate_impact_results.extend(moderate_impact)
        low_impact_results.extend(low_impact)
        modifier_impact_results.extend(modifier_impact)
        chromosome_data[chrom].append(combined_result)

results.sort(key=lambda x: (x[0], x[1]))
high_impact_results.sort(key=lambda x: (x[0], x[1]))
moderate_impact_results.sort(key=lambda x: (x[0], x[1]))
low_impact_results.sort(key=lambda x: (x[0], x[1]))
modifier_impact_results.sort(key=lambda x: (x[0], x[1]))

header_all = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
              "Chi2", "P-value", "CSQ", "Info", "Test Type"]
header_impact = ["Chromosome", "Position", "Ref", "Alt", "AC1", "AN1", "AC2", "AN2",
                 "Chi2", "P-value", "Info", "Test Type", "Gene Name", "Variant Type", "Impact"]

with open(output_csv, "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(header_all)
    csvwriter.writerows(results)
print(f"Results have been saved to {output_csv}")

for path, data, label in [
        (high_impact_csv, high_impact_results, "High"),
        (moderate_impact_csv, moderate_impact_results, "Moderate"),
        (low_impact_csv, low_impact_results, "Low"),
        (modifier_impact_csv, modifier_impact_results, "Modifier")]:
    with open(path, "w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(header_impact)
        csvwriter.writerows(data)
    print(f"{label} impact results have been saved to {path}")

for chrom in required_chromosomes:
    if chrom in chromosome_data:
        data = chromosome_data[chrom]
        data.sort(key=lambda x: x[1])
        chrom_csv = output_dir + f"NSCP_vs_SCAP-YAO-vep-output_{chrom}.csv"
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
