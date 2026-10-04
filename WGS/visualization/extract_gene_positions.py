"""
extract_gene_positions.py
Look up the coordinates (chromosome, start, end) of a list of genes in the
GFF3 annotations of the three reference genomes, and write a combined CSV:

    Gene | <ref1>-chr | <ref1>-start | <ref1>-end | <ref2>-... | <ref3>-...

Used to mark gene positions on chromosome ideograms across references.
"""

import pandas as pd

# --- User-configurable section ------------------------------------------------
gff_files = [
    "/path/to/anno/T2T-CHM13v2.0.sorted.gff3",
    "/path/to/anno/hg38.sorted.gff3",
    "/path/to/anno/T2T-Yao-hp.v1.1.sorted.gff3",
]
genome_labels = ["CHM13", "hg38", "YAO"]

# Genes to look up
genes_to_find = [
    "INTS11", "TTC22", "KLHDC8A", "IL15RA", "ZNF438", "LINC02666", "OR4D10",
    "VSIG10L2", "MAPKAPK5", "WDR95P", "TEX21P", "ZC2HC1C", "ATXN3", "IGF1R",
    "DCTPP1", "KIAA0753", "LGALS9B", "LINC02074", "ZNF527", "TMEM143",
    "LINC01722", "KAT14", "SNHG17", "IFNAR2", "SLC5A4", "MCM5", "CNTN4-AS1",
    "TTLL3", "TAMM41", "LINC00971", "LINC02614", "GRK4", "SERBP1P6", "CASC15",
    "FAM221A",
]
output_csv = "gene_positions.csv"
# -----------------------------------------------------------------------------


# Parse each GFF3 and record coordinates of genes in the target list
gene_info = {gene: {} for gene in genes_to_find}

for gff_file, label in zip(gff_files, genome_labels):
    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:            # malformed / truncated line — skip
                continue
            if parts[2] == "gene":
                info_fields = parts[8].split(";")
                gene_name = None
                for field in info_fields:
                    if field.startswith("Name="):
                        gene_name = field.split("=", 1)[1]
                if gene_name in genes_to_find:
                    gene_info[gene_name].update({
                        f"{label}-chr": parts[0],
                        f"{label}-start": parts[3],
                        f"{label}-end": parts[4],
                    })

# Build the wide-format table
data = []
for gene, info in gene_info.items():
    row = [gene]
    for label in genome_labels:
        row.extend([
            info.get(f"{label}-chr", "Not found"),
            info.get(f"{label}-start", "Not found"),
            info.get(f"{label}-end", "Not found"),
        ])
    data.append(row)

columns = ["Gene"]
for label in genome_labels:
    columns.extend([f"{label}-chr", f"{label}-start", f"{label}-end"])

df = pd.DataFrame(data, columns=columns)
df.to_csv(output_csv, index=False)
print(f"Gene positions written to {output_csv}")

# NOTE: for hg38 the "Name=" attribute of some GFF3 versions differs — genes
# not found for a reference show as "Not found"; adjust the attribute parsing
# if your annotation uses different keys.
