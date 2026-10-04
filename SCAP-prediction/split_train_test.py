"""
split_train_test.py

Split the LASSO-screened dataset into a training set and a held-out test set
with stratified random sampling, keeping the SCAP positive/negative ratio
identical in both partitions (80% / 20% by default).

The test set is NEVER used for any feature screening or model tuning in the
downstream pipeline — it only enters the final performance evaluation, and
its feature columns are aligned to whatever subset the training pipeline
selects.
"""

import pandas as pd
import os
from sklearn.model_selection import train_test_split

# --- User-configurable section ------------------------------------------------
# Input: the LASSO-reduced dataset exported by lasso_selection_glmnet.R
input_filepath = "/path/to/scap_genotype/Final_Data_For_Python_XGBoost.csv"
# Output directory for the split files
output_directory = "/path/to/scap_genotype/"
# Test fraction (0.2 = 80/20 split; use 0.3 for a 70/30 split)
test_ratio = 0.2
# Reproducibility
random_seed = 42
# Identifier / label columns
id_col = "ID"
label_col = "SCAP"
# ------------------------------------------------------------------------------


def split_clinical_dataset(input_file, output_dir, test_ratio=0.2,
                           random_seed=42):
    """Stratified train/test split preserving the class ratio."""
    print(f"1. Loading dataset: {input_file}")
    if not os.path.exists(input_file):
        print("File not found, please check the path.")
        return
    df = pd.read_csv(input_file)

    scap_counts = df[label_col].value_counts()
    print(f"\n[Dataset overview]")
    print(f"Total samples: {df.shape[0]}")
    print(f"Non-severe ({label_col}=0): {scap_counts.get(0, 0)}")
    print(f"Severe     ({label_col}=1): {scap_counts.get(1, 0)}")

    print(f"\n2. Performing stratified {int((1 - test_ratio) * 100)}/"
          f"{int(test_ratio * 100)} split...")
    train_df, test_df = train_test_split(
        df,
        test_size=test_ratio,
        random_state=random_seed,
        stratify=df[label_col]
    )

    print(f"\n[Split result]")
    print(f"--> Training set: {train_df.shape[0]} samples "
          f"({int(train_df[label_col].sum())} severe / "
          f"{train_df.shape[0] - int(train_df[label_col].sum())} non-severe)")
    print(f"--> Test set:     {test_df.shape[0]} samples "
          f"({int(test_df[label_col].sum())} severe / "
          f"{test_df.shape[0] - int(test_df[label_col].sum())} non-severe)")

    print("\n3. Saving split datasets...")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    train_output_path = os.path.join(output_dir, "Train_Data_80.csv")
    test_output_path = os.path.join(output_dir, "Test_Data_20.csv")
    train_df.to_csv(train_output_path, index=False)
    test_df.to_csv(test_output_path, index=False)
    print(f"Training set saved to: {train_output_path}")
    print(f"Test set saved to:     {test_output_path}")


if __name__ == "__main__":
    split_clinical_dataset(input_filepath, output_directory,
                           test_ratio=test_ratio, random_seed=random_seed)
