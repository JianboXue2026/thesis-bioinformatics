"""
consensus_feature_selection.py

Multi-model consensus feature screening on the training set only, after the
LASSO pre-filter. Three algorithms with different inductive biases rank the
candidate features in parallel (threads):

  Model A  Random Forest        - Gini-impurity importance, captures
                                  high-order non-linear interactions
  Model B  SVM-RFE              - maximum-margin hyperplane, robust in small
                                  samples, eliminates weakest feature first
  Model C  XGBoost + SHAP       - gradient-boosted trees with game-theoretic
                                  marginal contributions (mean |SHAP|)

Results per model are ranked; the Top-K sets are visualised with a Venn
diagram and combined into a final panel by one of two strategies:

  selection_mode = "pairwise"    (main pipeline) — save the three pairwise
      intersections, the strict 3-model intersection and the union of all
      pairwise intersections (features selected by >= 2 models). The final
      model is trained on Train/Test_Subset_Union_of_Pairs.csv.
  selection_mode = "topk_union"  (sensitivity analysis) — save the union of
      the three Top-K lists (feature-count sensitivity: top 20/25/30/35...).

The TEST SET is never used for ranking — it is only re-shaped to the selected
feature columns afterwards so both partitions stay aligned.

Outputs (per run, in <output_base_dir>/Consensus_Results_<timestamp>/):
  Result_*_Importance/Ranking.csv   full per-model rankings
  Result_Top<k>_*_Features.csv      (topk_union mode) per-model Top-K lists
  FigA_SHAP_Summary_Plot.pdf        SHAP beeswarm of the screening model
  FigB_Venn_Diagram_Intersection.pdf
  Train/Test_Subset_*.csv           aligned train/test subsets per set
"""

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
import shap
import concurrent.futures
import time
import datetime
import os

# --- User-configurable section ------------------------------------------------
# Training set (feature screening uses ONLY this partition)
train_path = "/path/to/scap_genotype/Train_Data_80.csv"
# Held-out test set (aligned extraction only, never screened)
test_path = "/path/to/scap_genotype/Test_Data_20.csv"
# Base output directory (a timestamped subdirectory is created per run)
output_base_dir = "/path/to/scap_genotype/"
# Number of top features per model entering the Venn diagram / sets
top_k = 20
# Consensus strategy: "pairwise" (main pipeline) or "topk_union" (sensitivity)
selection_mode = "pairwise"
# Identifier / label columns of the input table
id_col = "ID"
label_col = "SCAP"
# ------------------------------------------------------------------------------

# Publication-style vector figure settings
plt.rcParams["font.sans-serif"] = ["SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def run_rf(X, y, features):
    """Thread 1: Random Forest feature importance."""
    print("[Thread 1] Random Forest importance...")
    rf = RandomForestClassifier(n_estimators=300, random_state=42,
                                class_weight="balanced", n_jobs=-1)
    rf.fit(X, y)
    df_rf = pd.DataFrame(
        {"Feature": features, "RF_Score": rf.feature_importances_}
    ).sort_values(by="RF_Score", ascending=False).reset_index(drop=True)
    df_rf["RF_Rank"] = df_rf.index + 1
    return df_rf


def run_svm_rfe(X_scaled, y, features):
    """Thread 2: SVM recursive feature elimination (linear kernel)."""
    print("[Thread 2] SVM-RFE ranking...")
    svc = SVC(kernel="linear", random_state=42, class_weight="balanced")
    rfe = RFE(estimator=svc, n_features_to_select=1, step=1)
    rfe.fit(X_scaled, y)
    df_svm = pd.DataFrame(
        {"Feature": features, "SVM_Rank": rfe.ranking_}
    ).sort_values(by="SVM_Rank", ascending=True).reset_index(drop=True)
    return df_svm


def run_xgboost_shap(X, y, features):
    """Thread 3: XGBoost fit + mean |SHAP| importance."""
    print("[Thread 3] XGBoost + SHAP importance...")
    xgb_model = xgb.XGBClassifier(objective="binary:logistic",
                                  eval_metric="logloss", max_depth=4,
                                  learning_rate=0.05, n_estimators=200,
                                  random_state=42, n_jobs=-1)
    xgb_model.fit(X, y)
    explainer = shap.TreeExplainer(xgb_model)
    shap_values = explainer.shap_values(X)
    shap_importance = np.abs(shap_values).mean(axis=0)
    df_xgb = pd.DataFrame(
        {"Feature": features, "XGB_SHAP_Score": shap_importance}
    ).sort_values(by="XGB_SHAP_Score", ascending=False).reset_index(drop=True)
    df_xgb["XGB_Rank"] = df_xgb.index + 1
    return df_xgb, explainer, shap_values


def extract_and_save_subset(source_df, feature_set, output_path,
                            dataset_name, set_name):
    """Extract ID + label + feature_set columns from a dataframe and save."""
    if len(feature_set) == 0:
        print(f"  [Warning] {set_name} is empty, skipped for {dataset_name}.")
        return
    target_cols = [id_col, label_col] + list(feature_set)
    # Defensive: keep only columns that actually exist
    valid_cols = [col for col in target_cols if col in source_df.columns]
    subset_df = source_df[valid_cols]
    subset_df.to_csv(output_path, index=False)
    print(f"  Extracted [{set_name}] -> {os.path.basename(output_path)} "
          f"({len(valid_cols) - 2} features) for {dataset_name}")


def main(train_path, test_path, output_base_dir, top_k=20,
         selection_mode="pairwise"):
    start_time = time.time()
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(output_base_dir, f"Consensus_Results_{timestamp}")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("--- Stage 1: loading the training set ---")
    df_train = pd.read_csv(train_path)
    X_train = df_train.drop(columns=[id_col, label_col])
    y_train = df_train[label_col]
    features = X_train.columns.tolist()
    # SVM-RFE needs standardised features
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X_train), columns=features)
    print(f"Training set loaded: {X_train.shape[0]} samples, "
          f"{X_train.shape[1]} candidate features.")

    print("\n--- Stage 2: parallel multi-model screening (training set only) ---")
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        future_rf = executor.submit(run_rf, X_train, y_train, features)
        future_svm = executor.submit(run_svm_rfe, X_scaled, y_train, features)
        future_xgb = executor.submit(run_xgboost_shap, X_train, y_train,
                                     features)
        df_rf = future_rf.result()
        df_svm = future_svm.result()
        df_xgb, explainer, shap_values = future_xgb.result()

    print(f"\n--- Stage 3: saving per-model rankings and the SHAP figure ---")
    # Full rankings (kept for reference)
    df_rf.to_csv(os.path.join(output_dir, "Result_1_RF_Importance.csv"),
                 index=False)
    df_svm.to_csv(os.path.join(output_dir, "Result_2_SVM_RFE_Ranking.csv"),
                  index=False)
    df_xgb.to_csv(os.path.join(output_dir, "Result_3_XGBoost_SHAP_Ranking.csv"),
                  index=False)
    if selection_mode == "topk_union":
        df_rf.head(top_k).to_csv(
            os.path.join(output_dir, f"Result_Top{top_k}_RF_Features.csv"),
            index=False)
        df_svm.head(top_k).to_csv(
            os.path.join(output_dir, f"Result_Top{top_k}_SVM_Features.csv"),
            index=False)
        df_xgb.head(top_k).to_csv(
            os.path.join(output_dir, f"Result_Top{top_k}_XGBoost_Features.csv"),
            index=False)
    plt.figure(figsize=(10, 8))
    shap.summary_plot(shap_values, X_train, show=False, max_display=15)
    plt.title("XGBoost SHAP Summary Plot (Top 15 Features in Training Set)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "FigA_SHAP_Summary_Plot.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()

    print(f"\n--- Stage 4: Top {top_k} sets, Venn diagram and consensus sets ---")
    set_rf = set(df_rf["Feature"].head(top_k))
    set_svm = set(df_svm["Feature"].head(top_k))
    set_xgb = set(df_xgb["Feature"].head(top_k))

    plt.figure(figsize=(8, 8))
    venn3(subsets=(set_rf, set_svm, set_xgb),
          set_labels=("Random Forest", "SVM-RFE", "XGBoost (SHAP)"),
          set_colors=("#FF9999", "#99CCFF", "#99FF99"), alpha=0.7)
    plt.title(f"Intersection of Top {top_k} Features Across 3 Models",
              fontsize=16)
    plt.savefig(os.path.join(output_dir, "FigB_Venn_Diagram_Intersection.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()

    # Build the consensus feature sets to export
    if selection_mode == "pairwise":
        intersect_rf_svm = set_rf.intersection(set_svm)
        intersect_rf_xgb = set_rf.intersection(set_xgb)
        intersect_svm_xgb = set_svm.intersection(set_xgb)
        final_intersection = intersect_rf_svm.intersection(set_xgb)
        # Union of pairwise intersections = features chosen by >= 2 models
        union_of_pairs = (intersect_rf_svm.union(intersect_rf_xgb)
                          .union(intersect_svm_xgb))
        export_sets = [
            ("Pair1_RF_SVM", intersect_rf_svm),
            ("Pair2_RF_XGB", intersect_rf_xgb),
            ("Pair3_SVM_XGB", intersect_svm_xgb),
            ("Intersection_All3", final_intersection),
            ("Union_of_Pairs", union_of_pairs),
        ]
    elif selection_mode == "topk_union":
        topk_union = set_rf.union(set_svm).union(set_xgb)
        export_sets = [("Top" + str(top_k) + "_Union", topk_union)]
    else:
        raise ValueError(f"Unknown selection_mode: {selection_mode}")

    print(f"\n--- Stage 5: exporting train subsets per consensus set ---")
    for name, fset in export_sets:
        extract_and_save_subset(
            df_train, fset,
            os.path.join(output_dir, f"Train_Subset_{name}.csv"),
            "training set", name)

    print(f"\n--- Stage 6: aligning the test set to the same features ---")
    if os.path.exists(test_path):
        df_test = pd.read_csv(test_path)
        for name, fset in export_sets:
            extract_and_save_subset(
                df_test, fset,
                os.path.join(output_dir, f"Test_Subset_{name}.csv"),
                "test set", name)
        print("Test set alignment done.")
    else:
        print(f"[Warning] Test set {test_path} not found, skipped.")

    print("\n=============================================")
    print("[Summary]")
    print(f"- Top-K threshold: {top_k}")
    print(f"- Selection mode:  {selection_mode}")
    for name, fset in export_sets:
        print(f"- {name}: {len(fset)} features")
    print("=============================================")
    elapsed_time = time.time() - start_time
    print(f"\nAll done in {elapsed_time:.2f} s. Results in: {output_dir}")


if __name__ == "__main__":
    main(train_path=train_path, test_path=test_path,
         output_base_dir=output_base_dir, top_k=top_k,
         selection_mode=selection_mode)
