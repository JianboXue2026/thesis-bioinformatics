"""
plot_confusion_corr.py

Supplementary evaluation plots from the saved model assets and the cleaned
test subset (both produced by train_final_model.py):

  1. Confusion matrix on the test set (prediction at the model's default
     0.5 cut-off) rendered as a heatmap.
  2. Feature correlation heatmap of the panel (top-20 features by XGBoost
     importance when available, otherwise the first 20 columns) — visual
     check of residual collinearity in the final panel.
"""

import pandas as pd
import numpy as np
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

# --- User-configurable section ------------------------------------------------
# Model assets and cleaned test subset (both from the same run!)
# NOTE: Test_Subset_Cleaned_Final.csv is only written when remove_collinear =
# True in train_final_model.py. For the main pipeline (remove_collinear =
# False) point data_path to Test_Subset_Union_of_Pairs.csv instead.
model_path = "/path/to/scap_genotype/Consensus_Results_YYYYMMDD_HHMMSS/Final_SCAP_Model_Assets.pkl"
data_path = "/path/to/scap_genotype/Consensus_Results_YYYYMMDD_HHMMSS/Evaluation/Test_Subset_Cleaned_Final.csv"
# Output directory for the two PDFs
output_dir = "/path/to/scap_genotype/Consensus_Results_YYYYMMDD_HHMMSS/Evaluation/"
# Column names
id_column = "ID"
target_column = "SCAP"
# ------------------------------------------------------------------------------


def main():
    # 1. Load the test data and the model assets
    print(f"Loading data from {data_path}...")
    df = pd.read_csv(data_path)
    print(f"Loading model from {model_path}...")
    assets = joblib.load(model_path)

    # 2. Locate a usable model object inside the assets dict
    # Prefer the XGBoost model explicitly: the assets dict also contains the
    # sklearn logistic regression (trained on STANDARDISED features) — picking
    # it by accident and calling predict on RAW features would silently
    # produce wrong confusion matrices.
    model = None
    feature_names = None
    if isinstance(assets, dict):
        for f_key in ["selected_features", "feature_names", "features",
                      "genes"]:
            if f_key in assets:
                feature_names = assets[f_key]
                break
        if "xgb_model" in assets:
            model = assets["xgb_model"]
            print("Model object found (key: 'xgb_model').")
        else:
            for k, v in assets.items():
                if hasattr(v, "predict"):
                    model = v
                    print(f"Model object found (key: '{k}').")
                    break
            if model is None:
                raise ValueError("No object with a predict method found in assets!")
    else:
        model = assets
        if not hasattr(model, "predict"):
            raise ValueError("Loaded object has no predict method!")

    # 3. Feature matrix / labels (align column order to the stored list)
    X_test = df.drop(columns=[id_column, target_column])
    y_true = df[target_column]
    if feature_names is not None:
        X_test = X_test[feature_names]

    print("Predicting...")
    y_pred = model.predict(X_test)

    # 4. Plot 1: confusion matrix
    print("Saving confusion matrix (PDF)...")
    plt.figure(figsize=(8, 6))
    cm = confusion_matrix(y_true, y_pred)
    sns.heatmap(cm, annot=True, fmt="d", cmap="Blues",
                xticklabels=["Non-SCAP (0)", "SCAP (1)"],
                yticklabels=["Non-SCAP (0)", "SCAP (1)"])
    plt.title("Confusion Matrix: Predicted vs Actual", fontsize=14, pad=15)
    plt.xlabel("Predicted Label", fontsize=12)
    plt.ylabel("True Label", fontsize=12)
    cm_output = output_dir + "Confusion_Matrix.pdf"
    plt.savefig(cm_output, format="pdf", bbox_inches="tight")
    plt.close()

    # 5. Plot 2: feature correlation heatmap (top-20 by importance if the
    #    model exposes feature_importances_)
    print("Saving feature correlation heatmap (PDF)...")
    plt.figure(figsize=(10, 8))
    try:
        if hasattr(model, "feature_importances_"):
            importances = model.feature_importances_
            indices = np.argsort(importances)[-20:][::-1]
            top_features = X_test.columns[indices].tolist()
            corr_data = X_test[top_features]
        else:
            corr_data = X_test.iloc[:, :20]
    except Exception:
        corr_data = X_test.iloc[:, :20]
    corr_matrix = corr_data.corr()
    mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
    sns.heatmap(corr_matrix, mask=mask, annot=True, cmap="RdBu_r", center=0,
                fmt=".2f", annot_kws={"size": 8})
    plt.title("Correlation Heatmap", fontsize=14, pad=15)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    corr_output = output_dir + "Feature_Correlation_Heatmap.pdf"
    plt.savefig(corr_output, format="pdf", bbox_inches="tight")
    plt.close()

    print(f"\nDone. Two PDFs written to {output_dir}:")
    print(f"  1. {cm_output}")
    print(f"  2. {corr_output}")


if __name__ == "__main__":
    main()
