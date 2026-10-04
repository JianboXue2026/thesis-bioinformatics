"""
predict_new_samples.py

Single / batch inference with the trained SCAP panel assets
(Final_SCAP_Model_Assets.pkl produced by train_final_model.py).

IMPORTANT — feature scaling convention:
  * XGBoost was trained on the RAW 0/1 features       -> no transformation
  * Logistic regression was trained on STANDARDISED
    features (the StandardScaler fitted on the training set is stored in
    the assets)                                        -> scaler.transform()

Input CSV must contain: an ID column plus one column per panel locus
(either the original 0/1 one-hot indicators, or the raw -1/0/1/2 genotype
coding if the panel was trained on that representation — the column NAMES
must match the stored gene list exactly).
"""

import pandas as pd
import joblib
import os

# --- User-configurable section ------------------------------------------------
# Model assets exported by train_final_model.py
assets_path = "/path/to/scap_genotype/Consensus_Results_YYYYMMDD_HHMMSS/Final_SCAP_Model_Assets.pkl"
# New patient data (CSV)
new_data = "/path/to/new_patients.csv"
# Output file for the prediction results
output_path = "/path/to/scap_genotype/Inference_Results.csv"
# Identifier column in the new-patient CSV
id_col = "ID"
# Decision threshold for the risk label (Youden cut-off from training;
# adjust to the value reported by train_final_model.py)
risk_threshold = 0.5
# ------------------------------------------------------------------------------


def run_inference(assets_path, new_patient_csv, output_path,
                  risk_threshold=0.5):
    # 1. Load the model assets
    print(f"Loading model assets: {assets_path}...")
    assets = joblib.load(assets_path)
    xgb_model = assets["xgb_model"]
    lr_model = assets["lr_model"]
    scaler = assets["scaler"]
    required_genes = assets["genes"]

    # 2. Read the new patient data and validate the feature columns
    df_new = pd.read_csv(new_patient_csv)
    missing_genes = [g for g in required_genes if g not in df_new.columns]
    if missing_genes:
        raise ValueError(
            f"Input CSV is missing required loci: {missing_genes}")
    X_new = df_new[required_genes]

    # 3. Predict — see the scaling note in the module docstring
    print(f"Computing SCAP probabilities for {len(df_new)} samples...")
    xgb_probs = xgb_model.predict_proba(X_new)[:, 1]
    X_new_scaled = scaler.transform(X_new)
    lr_probs = lr_model.predict_proba(X_new_scaled)[:, 1]

    # 4. Assemble and save the results
    results = pd.DataFrame({
        "Patient_ID": df_new[id_col],
        "XGBoost_Probability": xgb_probs,
        "Logistic_Probability": lr_probs,
        "Risk_Level": ["High Risk" if p > risk_threshold else "Low Risk"
                       for p in xgb_probs],
    })
    results.to_csv(output_path, index=False)
    print(f"Done. Predictions saved to: {output_path}")
    print(results.head())


if __name__ == "__main__":
    run_inference(assets_path, new_data, output_path,
                  risk_threshold=risk_threshold)
