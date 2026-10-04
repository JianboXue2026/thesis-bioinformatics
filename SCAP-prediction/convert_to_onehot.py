"""
convert_to_onehot.py

Convert a genotype table (one row per patient) into one-hot encoded features.

Input CSV layout:
    ID     : patient identifier
    SCAP   : phenotype label (0 = non-severe, 1 = severe CAP)
    YAO001, YAO002, ... : genotype columns coded as
        0 = wild type, 1 = heterozygous mutation, 2 = homozygous mutation,
       -1 = unknown / not called

Each genotype column is expanded into up to 4 binary indicator columns
(e.g. YAO001_0, YAO001_1, YAO001_2, YAO001_-1), so downstream penalised
models see mutually exclusive 0/1 dummy variables instead of an ordinal
coding.
"""

import pandas as pd
import os

# --- User-configurable section ------------------------------------------------
input_file = "/path/to/scap_genotype/SCAP_GenoType.csv"
output_file = "/path/to/scap_genotype/SCAP_GenoType_OneHot.csv"
# Metadata columns that must be carried through unchanged
meta_cols = ["ID", "SCAP"]
# ------------------------------------------------------------------------------


def convert_to_onehot(input_path, output_path):
    """Convert genotype data coded as -1/0/1/2 into one-hot (dummy) encoding."""
    print(f"1. Reading input file: {input_path}")
    if not os.path.exists(input_path):
        print("Input file not found, please check the path.")
        return
    df = pd.read_csv(input_path)

    # Separate metadata (ID, phenotype label) from genotype features
    meta_data = df[meta_cols]
    gene_data = df.drop(columns=meta_cols)
    print(f"Original feature dimension: {gene_data.shape[1]} genotype loci")

    # Cast to string so 0/1/2/-1 are treated as 4 independent categories
    # (prevents pandas from interpreting them as ordinal numbers)
    gene_data_str = gene_data.astype(str)

    print("\n2. Running one-hot encoding...")
    # dtype=int guarantees a numeric 0/1 matrix instead of True/False
    gene_data_onehot = pd.get_dummies(gene_data_str, dtype=int)
    print(f"Dimension after one-hot encoding: {gene_data_onehot.shape[1]} columns")

    print("\n3. Re-assembling metadata and encoded features...")
    final_df = pd.concat([meta_data, gene_data_onehot], axis=1)
    final_df.to_csv(output_path, index=False)
    print(f"Done. One-hot encoded data saved to: {output_path}")

    print("\n[Preview]")
    print(final_df.iloc[:, :10].head())


if __name__ == "__main__":
    convert_to_onehot(input_file, output_file)
