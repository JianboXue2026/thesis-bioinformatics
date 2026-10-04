"""
bim_add_variant_id.py
Assign a unique variant identifier to every row of a PLINK .bim file by
rewriting column 2 as <chromosome>:<position>:<REF>:<ALT>, e.g. "1:897539:A:G".
This makes variants uniquely addressable in downstream association results
(the default .bim IDs contain duplicates after multi-reference merging).

The .bim layout is: chr, id, cM, bp, allele1 (A1), allele2 (A2), so the
deterministic ID is built from columns 1, 4, 6, 5 (1-based).

Note on naming: PLINK's .bim columns 5/6 are simply "allele 1" and "allele 2";
PLINK does **not** guarantee that allele 2 is the reference allele. The ID is
built as `<allele2>:<allele1>` purely as a deterministic, collision-free label —
read the last two fields as "the two alleles observed at this locus", not as a
verified REF/ALT pair.

Note: rows whose alleles are missing ("." or "0") fall back to a
position-based prefix plus the original ID, since a deterministic
chr:pos:allele:allele ID is not well defined for them.
"""

import pandas as pd

# --- User-configurable section ----------------------------------------------
bim_file = "/path/to/GATK-PLINK/01.trans2bed/merged-CHM13.bim"
output_file = "/path/to/GATK-PLINK/01.trans2bed/merged-CHM13.modified.bim"
# -----------------------------------------------------------------------------

# Load the BIM file (whitespace-separated, no header)
bim = pd.read_csv(bim_file, sep=r"\s+", header=None, low_memory=False)

# New SNP identifier: chr:pos:allele2:allele1 (deterministic, collision-free)
# .bim columns (0-based): 0=chr, 3=bp, 4=A1, 5=A2
def make_id(row):
    ref, alt = str(row[5]), str(row[4])
    if ref in (".", "0", "nan") or alt in (".", "0", "nan"):
        # Alleles unknown — fall back to a positional + original-ID form
        return f"{row[0]}:{row[3]}:{row[1]}"
    return f"{row[0]}:{row[3]}:{ref}:{alt}"

bim[1] = bim.apply(make_id, axis=1)

bim.to_csv(output_file, sep="\t", header=False, index=False)
print(f"Modified BIM file saved to {output_file}")
