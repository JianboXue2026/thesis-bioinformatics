"""
train_final_model.py

Build and evaluate the final SCAP genotype-panel prediction model on the
consensus-screened feature subsets:

  - XGBoost                      : non-linear model, maximises discrimination
                                   (AUC); trained on the RAW 0/1 features
  - Multivariable logistic reg.  : linear model for clinical interpretability;
                                   trained on STANDARDISED features (sklearn)
                                   and on raw features via statsmodels (the
                                   statsmodels fit provides the coefficients
                                   and p-values used by the nomogram)

Evaluation on the held-out test set (Youden-optimal cut-off):
  AUC, Brier score, sensitivity, specificity, PPV, NPV.

Outputs (in <timestamped_results_dir>/Evaluation/):
  Evaluation_Fig1_ROC_Curve.pdf               ROC curves, both models
  Evaluation_Fig2_Calibration_Curve.pdf       calibration curves, both models
  Evaluation_Fig3_Decision_Curve_Analysis.pdf DCA (net benefit, computed
                                              manually - no external DCA lib)
  Evaluation_Fig4_Clinical_Nomogram.pdf       nomogram (points per locus ->
                                              total points -> SCAP probability)
  Clinical_Score_System_Text_Version.txt      text version of the score system
  Evaluation_Fig5_SHAP_Beeswarm.pdf           SHAP interpretation of XGBoost
  Train/Test_Subset_Cleaned_Final.csv         (if remove_collinear = True)
  ../Final_SCAP_Model_Assets.pkl              all fitted objects for inference
"""

# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_curve, auc, brier_score_loss,
                             confusion_matrix, classification_report)
from sklearn.calibration import calibration_curve
import xgboost as xgb
import os
import joblib
import shap
import statsmodels.api as sm

# --- User-configurable section ------------------------------------------------
# Directory created by consensus_feature_selection.py
timestamped_results_dir = "/path/to/scap_genotype/Consensus_Results_YYYYMMDD_HHMMSS"
# Feature-subset files inside that directory
# main pipeline: Train/Test_Subset_Union_of_Pairs.csv
# topk-union mode: Train/Test_Subset_Top20_Union.csv
train_subset_name = "Train_Subset_Union_of_Pairs.csv"
test_subset_name = "Test_Subset_Union_of_Pairs.csv"
# Data-cleaning switch: remove zero-variance features and features with
# pairwise |r| > 0.95 before modelling (used for the all-LASSO-features
# sensitivity analysis; keep False for the union-of-pairs main panel)
remove_collinear = False
corr_threshold = 0.95
# Column names
id_col = "ID"
label_col = "SCAP"
# ------------------------------------------------------------------------------

# Publication-style vector figure settings
plt.rcParams["font.sans-serif"] = ["SimHei", "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


def train_evaluate_final_panel(timestamped_results_dir):
    """Train the final XGBoost + logistic-regression panel models, evaluate
    them on the held-out test set, and export all figures + model assets."""
    print("=" * 58)
    print("Starting final panel training & evaluation...")
    print(f"Data source directory: {timestamped_results_dir}")
    print("=" * 58 + "\n")

    evaluation_dir = os.path.join(timestamped_results_dir, "Evaluation")
    if not os.path.exists(evaluation_dir):
        os.makedirs(evaluation_dir)
    train_subset_path = os.path.join(timestamped_results_dir,
                                     train_subset_name)
    test_subset_path = os.path.join(timestamped_results_dir,
                                    test_subset_name)

    # ------------------------------------------------------------------
    # Stage 1: load the screened data subsets
    # ------------------------------------------------------------------
    print("1. Loading screened data subsets...")
    if not os.path.exists(train_subset_path) or not os.path.exists(test_subset_path):
        print("Error: subset files not found, check the previous step.")
        return
    df_train = pd.read_csv(train_subset_path)
    df_test = pd.read_csv(test_subset_path)
    final_genes = [c for c in df_train.columns if c not in [id_col, label_col]]
    print(f"   Features entering the panel: {len(final_genes)}")
    X_train = df_train[final_genes]
    y_train = df_train[label_col]
    X_test = df_test[final_genes]
    y_test = df_test[label_col]

    # ------------------------------------------------------------------
    # Stage 1.5 (optional): wash-out of zero-variance and highly
    # collinear features (prevents a singular design matrix)
    # ------------------------------------------------------------------
    if remove_collinear:
        print("\n1.5 Data wash-out (zero variance / collinearity)...")
        variances = X_train.var()
        zero_var_cols = variances[variances == 0].index.tolist()
        if zero_var_cols:
            print(f"   [Drop] zero-variance features ({len(zero_var_cols)}): "
                  f"{zero_var_cols}")
        else:
            print("   [Check] no zero-variance features.")
        X_train_temp = X_train.drop(columns=zero_var_cols)
        corr_matrix = X_train_temp.corr().abs()
        upper = corr_matrix.where(
            np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        high_corr_cols = [c for c in upper.columns
                          if any(upper[c] > corr_threshold)]
        if high_corr_cols:
            print(f"   [Drop] collinear features |r| > {corr_threshold} "
                  f"({len(high_corr_cols)}): {high_corr_cols}")
        else:
            print(f"   [Check] no features with |r| > {corr_threshold}.")
        cols_to_drop = list(set(zero_var_cols + high_corr_cols))
        if cols_to_drop:
            X_train = X_train.drop(columns=cols_to_drop)
            X_test = X_test.drop(columns=cols_to_drop)
            final_genes = X_train.columns.tolist()
            print(f"   -> {len(final_genes)} independent features remain.")
        # Persist the cleaned datasets (ID + label + final features)
        df_train[[id_col, label_col] + final_genes].to_csv(
            os.path.join(evaluation_dir, "Train_Subset_Cleaned_Final.csv"),
            index=False)
        df_test[[id_col, label_col] + final_genes].to_csv(
            os.path.join(evaluation_dir, "Test_Subset_Cleaned_Final.csv"),
            index=False)
        print("   [Saved] cleaned train/test subsets -> Evaluation/")

    # Standardised copy for the (sklearn) logistic regression
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # ------------------------------------------------------------------
    # Stage 2: XGBoost (maximise discrimination)
    # ------------------------------------------------------------------
    print("\n2. Fitting XGBoost (discrimination-oriented)...")
    xgb_model = xgb.XGBClassifier(
        objective="binary:logistic", eval_metric="auc",
        max_depth=4, learning_rate=0.03, n_estimators=300,
        subsample=0.8, colsample_bytree=0.8,
        random_state=42, n_jobs=-1
    )
    xgb_model.fit(X_train, y_train)

    # ------------------------------------------------------------------
    # Stage 3: multivariable logistic regression (nomogram-oriented)
    # ------------------------------------------------------------------
    print("3. Fitting multivariable logistic regression "
          "(interpretability / nomogram)...")
    # statsmodels fit on raw features: coefficients and p-values feed the
    # nomogram; explicit constant term (intercept)
    X_train_logit_sm = sm.add_constant(X_train)
    X_test_logit_sm = sm.add_constant(X_test)
    try:
        logit_model_sm = sm.Logit(y_train, X_train_logit_sm).fit(
            method="bfgs", maxiter=1000, disp=0)
    except Exception:
        print("   [Warning] BFGS failed, falling back to lbfgs...")
        logit_model_sm = sm.Logit(y_train, X_train_logit_sm).fit(
            method="lbfgs", maxiter=1000, disp=0)
    # sklearn fit on standardised features: fast, class-weighted predictions
    lr_model = LogisticRegression(solver="liblinear", random_state=42,
                                  class_weight="balanced")
    lr_model.fit(X_train_scaled, y_train)

    # ------------------------------------------------------------------
    # Stage 4: full evaluation on the held-out test set
    # ------------------------------------------------------------------
    print("\n4. Evaluating on the test set...")

    def evaluate_on_test(model, X_input, y_true, model_name, is_sklearn=True):
        if is_sklearn:
            y_pred_proba = model.predict_proba(X_input)[:, 1]
        else:
            y_pred_proba = model.predict(X_input)
        fpr, tpr, thresholds = roc_curve(y_true, y_pred_proba)
        roc_auc = auc(fpr, tpr)
        brier = brier_score_loss(y_true, y_pred_proba)
        # Youden index: maximise sensitivity + specificity - 1
        youden_index = tpr - fpr
        best_threshold = thresholds[np.argmax(youden_index)]
        # Re-classify at the optimal cut-off (default would be 0.5)
        y_pred_at_best = (y_pred_proba > best_threshold).astype(int)
        cm = confusion_matrix(y_true, y_pred_at_best)
        tn, fp, fn, tp = cm.ravel()
        sensitivity = tp / (tp + fn)
        specificity = tn / (tn + fp)
        ppv = tp / (tp + fp) if (tp + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        print("\n" + "=" * 45)
        print(f"Test-set performance: {model_name}\n")
        print(f"[Discrimination]  AUC:         {roc_auc:.4f}")
        print(f"[Calibration]    Brier score: {brier:.4f}")
        print("=" * 45)
        print(f"Youden-optimal cut-off: {best_threshold:.4f}\n")
        print(f"Sensitivity (TPR): {sensitivity * 100:.1f}%")
        print(f"Specificity (TNR): {specificity * 100:.1f}%")
        print(f"PPV:               {ppv * 100:.1f}%")
        print(f"NPV:               {npv * 100:.1f}%")
        print("=" * 45)
        return {
            "Model": model_name, "AUC": roc_auc, "Brier": brier,
            "Cut-off": best_threshold,
            "Sensitivity_Youden": sensitivity,
            "Specificity_Youden": specificity,
            "PPV": ppv, "NPV": npv,
            "fpr": fpr, "tpr": tpr, "y_pred_proba": y_pred_proba,
        }

    # XGBoost evaluated on raw features; LR on standardised features
    xgb_results = evaluate_on_test(xgb_model, X_test, y_test,
                                   "XGBoost Final Model")
    lr_results = evaluate_on_test(lr_model, X_test_scaled, y_test,
                                  "Logistic Regression (Scaled)")

    # ------------------------------------------------------------------
    # Stage 5: publication figures (vector PDF)
    # ------------------------------------------------------------------
    print("\n5. Generating vector figures (PDF)...")

    # Fig 1: combined ROC curve with the Youden cut-off marked
    plt.figure(figsize=(7, 6))
    plt.plot(xgb_results["fpr"], xgb_results["tpr"], color="#d62728", lw=2.5,
             label=f"XGBoost (AUC = {xgb_results['AUC']:.3f})")
    plt.plot(lr_results["fpr"], lr_results["tpr"], color="#1f77b4", lw=2,
             linestyle="--",
             label=f"Logistic Regression (AUC = {lr_results['AUC']:.3f})")
    plt.plot([0, 1], [0, 1], color="gray", lw=1, linestyle="-")
    cut_idx = np.argmax(xgb_results["tpr"] - xgb_results["fpr"])
    plt.scatter(xgb_results["fpr"][cut_idx], xgb_results["tpr"][cut_idx],
                color="darkred", marker="o", s=100, edgecolors="black",
                label="Youden Cut-off (XGB)")
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel("False Positive Rate", fontsize=12)
    plt.ylabel("True Positive Rate", fontsize=12)
    plt.title("Gene Panel Predictive Performance (ROC)", fontsize=14,
              fontweight="bold")
    plt.legend(loc="lower right", fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(evaluation_dir, "Evaluation_Fig1_ROC_Curve.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()

    # Fig 2: calibration curves (goodness of fit)
    print("   Calibration curves...")
    prob_true_lr, prob_pred_lr = calibration_curve(
        y_test, lr_results["y_pred_proba"], n_bins=10, strategy="uniform")
    prob_true_xgb, prob_pred_xgb = calibration_curve(
        y_test, xgb_results["y_pred_proba"], n_bins=10, strategy="uniform")
    plt.figure(figsize=(7, 6))
    plt.plot(prob_pred_lr, prob_true_lr, "s-", color="#1f77b4", marker="s",
             markersize=6,
             label=f"Logistic Regression (Brier={lr_results['Brier']:.3f})")
    plt.plot(prob_pred_xgb, prob_true_xgb, "o-", color="#d62728", marker="o",
             markersize=6,
             label=f"XGBoost (Brier={xgb_results['Brier']:.3f})")
    plt.plot([0, 1], [0, 1], "k:", label="Perfectly calibrated (Ideal)")
    plt.ylabel("Actual observed probability", fontsize=12)
    plt.xlabel("Model-predicted SCAP probability", fontsize=12)
    plt.title("Model Calibration Curve", fontsize=14, fontweight="bold")
    plt.legend(loc="lower right")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(evaluation_dir,
                             "Evaluation_Fig2_Calibration_Curve.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()

    # Fig 3: decision curve analysis — net benefit computed manually,
    # no external DCA package required
    print("   Decision curve analysis (DCA)...")

    def calculate_net_benefit(y_true, y_pred_proba, threshold):
        """Net benefit at a given threshold probability."""
        tn, fp, fn, tp = confusion_matrix(
            y_true, (y_pred_proba > threshold).astype(int)).ravel()
        n_total = tn + fp + fn + tp
        return (tp / n_total) - (fp / n_total) * (
            threshold / (1 - threshold + 1e-9))

    plt.figure(figsize=(7, 6))
    thresholds_dca = np.arange(0.01, 1.0, 0.01)
    nb_xgb, nb_lr = [], []
    for t in thresholds_dca:
        nb_xgb.append(calculate_net_benefit(y_test,
                                            xgb_results["y_pred_proba"], t))
        nb_lr.append(calculate_net_benefit(y_test,
                                           lr_results["y_pred_proba"], t))
    scap_pos_rate = y_test.sum() / len(y_test)
    # "Treat all" and "treat none" reference strategies
    nb_treat_all = [scap_pos_rate - (1 - scap_pos_rate) * (t / (1 - t))
                    for t in thresholds_dca]
    nb_treat_none = [0] * len(thresholds_dca)
    plt.plot(thresholds_dca, nb_xgb, color="#d62728", lw=2.5,
             label="Gene Panel (XGBoost)")
    plt.plot(thresholds_dca, nb_lr, color="#1f77b4", lw=2, linestyle="--",
             label="Gene Panel (Logistic Regression)")
    plt.plot(thresholds_dca, nb_treat_all, color="gray", lw=1, linestyle=":",
             label="Treat All patients")
    plt.plot(thresholds_dca, nb_treat_none, color="black", lw=1,
             label="Treat None patients")
    plt.xlim([0.0, 1.0])
    plt.ylim([-0.05, max(max(nb_xgb), max(nb_lr)) * 1.2])
    plt.xlabel("Threshold Probability", fontsize=12)
    plt.ylabel("Net Benefit", fontsize=12)
    plt.title("Clinical Decision Curve Analysis (DCA)", fontsize=14,
              fontweight="bold")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(
        evaluation_dir, "Evaluation_Fig3_Decision_Curve_Analysis.pdf"),
        format="pdf", bbox_inches="tight")
    plt.close()

    # Fig 4: clinical nomogram + text version of the score system,
    # derived from the statsmodels logistic-regression coefficients
    print("   Clinical nomogram and text score system...")
    logit_summary = logit_model_sm.summary2().tables[1]
    if "Coef." in logit_summary.columns:
        logit_summary.rename(columns={"Coef.": "Coef", "coef": "Coef"},
                             inplace=True)
    gene_coefs = logit_summary.drop("const").sort_values(by="Coef",
                                                        ascending=False)
    num_vars = len(gene_coefs)
    # +3 rows: points axis, total-points axis, probability axis
    fig = plt.figure(figsize=(12, max(9, num_vars * 1.0)))
    gs = gridspec.GridSpec(num_vars + 3, 1)

    # Top axis: 0-100 points
    ax_points = plt.subplot(gs[0])
    ax_points.plot([0, 100], [0, 0], "k-", lw=1.5)
    points_ticks = np.arange(0, 101, 10)
    ax_points.set_xticks(points_ticks)
    ax_points.set_xticklabels(points_ticks)
    ax_points.set_title("Points", fontsize=12, pad=10)
    ax_points.set_xlim([-5, 105])
    ax_points.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax_points.spines[sp].set_visible(False)
    for p in points_ticks:
        ax_points.axvline(p, ymin=0, ymax=0.5, color="k", lw=1)

    # Scale: the largest |coefficient| maps to 100 points
    max_coef_abs = gene_coefs["Coef"].abs().max()
    scaling_factor = 100 / max_coef_abs
    text_report = [
        "=" * 50,
        "SCAP targeted-genotype clinical score system",
        "=" * 50,
        "Instructions: sum the points of each locus genotype.\n",
        "[Part 1: points per locus]",
    ]

    # One axis row per locus
    for i, (feat_name, row) in enumerate(gene_coefs.iterrows()):
        ax_var = plt.subplot(gs[i + 1])
        coef_val = row["Coef"]
        mutation_points = coef_val * scaling_factor
        if coef_val > 0:
            p_low, p_high = 0, mutation_points
            label_low, label_high = "Wild (0)", "Mutated (1)"
            text_report.append(
                f"Locus {feat_name:10}: mutated = {round(mutation_points):3} "
                f"points | wild = 0 points")
        else:
            p_low, p_high = 0, np.abs(mutation_points)
            label_low, label_high = "Mutated (1)", "Wild (0)"
            text_report.append(
                f"Locus {feat_name:10}: wild = {round(np.abs(mutation_points)):3} "
                f"points | mutated = 0 points (protective)")
        ax_var.plot([p_low, p_high], [0, 0], "k-", lw=1.5)
        ax_var.axvline(p_low, ymin=0, ymax=0.5, color="k", lw=1)
        ax_var.axvline(p_high, ymin=0, ymax=0.5, color="k", lw=1)
        ax_var.set_ylabel(f"{feat_name}\n(p={row['P>|z|']:.3f})", rotation=0,
                          labelpad=30, ha="right", va="center", fontsize=9)
        ax_var.text(p_low, 0.1, label_low, ha="center", fontsize=8)
        ax_var.text(p_high, 0.1, label_high, ha="center", fontsize=8)
        ax_var.set_xlim([-5, 105])
        ax_var.set_yticks([])
        ax_var.set_xticks([])
        for sp in ["top", "right", "left", "bottom"]:
            ax_var.spines[sp].set_visible(False)

    # Total points axis
    max_total_points = np.abs(gene_coefs["Coef"]).sum() * scaling_factor
    points_range = np.arange(0, max_total_points + 50, 50)
    ax_total = plt.subplot(gs[num_vars + 1])
    ax_total.plot([0, max(points_range)], [0, 0], "k-", lw=2)
    ax_total.set_xticks(points_range)
    ax_total.set_title("Total Points", fontsize=12, pad=10)
    for sp in ["top", "right", "left"]:
        ax_total.spines[sp].set_visible(False)
    for p in points_range:
        ax_total.axvline(p, ymin=0, ymax=0.5, color="k", lw=1)
    ax_total.set_xlim([-5, points_range[-1] + 5])
    ax_total.set_yticks([])

    # Probability axis: total points -> SCAP probability via the log-odds
    # transform with the fitted intercept (must align with Total Points)
    intercept = logit_model_sm.params["const"]
    targets_prob = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.85, 0.95])
    ax_prob = plt.subplot(gs[num_vars + 2])
    ax_prob.plot([0, max(points_range)], [0, 0], "k-", lw=2)
    prob_ticks, prob_labels = [], []
    text_report.append("\n[Part 2: total points -> SCAP probability]")
    for p in targets_prob:
        log_odds = np.log(p / (1 - p))
        coef_sum_needed = log_odds - intercept
        total_points_for_p = coef_sum_needed * scaling_factor
        text_report.append(
            f"Total points {round(total_points_for_p):4}  ==>  "
            f"SCAP probability ~ {p * 100:2.0f}%")
        if -50 <= total_points_for_p <= max(points_range) + 50:
            prob_ticks.append(total_points_for_p)
            prob_labels.append(f"{p * 100:.0f}%")
            ax_prob.axvline(total_points_for_p, ymin=0, ymax=0.5, color="k",
                            lw=1)
    ax_prob.set_xticks(prob_ticks)
    ax_prob.set_xticklabels(prob_labels, fontsize=10, fontweight="bold")
    ax_prob.set_title("SCAP Probability", fontsize=12, pad=10)
    ax_prob.set_xlim([-5, points_range[-1] + 5])
    ax_prob.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax_prob.spines[sp].set_visible(False)
    plt.suptitle("Clinical Nomogram to Predict Severe Pneumonia (SCAP)",
                 fontsize=16, fontweight="bold", y=0.98)
    # Extra row padding keeps neighbouring tick labels from overlapping
    plt.tight_layout(rect=[0, 0.03, 1, 0.95], h_pad=1.5)
    plt.savefig(os.path.join(evaluation_dir,
                             "Evaluation_Fig4_Clinical_Nomogram.pdf"),
                format="pdf", bbox_inches="tight")
    plt.close()

    report_text = "\n".join(text_report)
    print("\n" + report_text + "\n")
    with open(os.path.join(evaluation_dir,
                           "Clinical_Score_System_Text_Version.txt"), "w",
              encoding="utf-8") as f:
        f.write(report_text)

    # Fig 5: SHAP beeswarm of the final XGBoost
    print("6. Generating the XGBoost SHAP beeswarm plot...")
    explainer = shap.TreeExplainer(xgb_model)
    shap_values = explainer(X_train)
    plt.figure(figsize=(10, 8))
    shap.plots.beeswarm(shap_values, max_display=len(final_genes), show=False)
    plt.title("SHAP Feature Interpretation (Impact on SCAP Risk)",
              fontsize=14, pad=20)
    plt.savefig(os.path.join(evaluation_dir,
                             "Evaluation_Fig5_SHAP_Beeswarm.pdf"),
                bbox_inches="tight")
    plt.close()

    # ------------------------------------------------------------------
    # Stage 6: persist model assets for inference
    # ------------------------------------------------------------------
    print("\n7. Saving model assets...")
    model_assets = {
        "xgb_model": xgb_model,
        "lr_model": lr_model,
        "scaler": scaler,
        "genes": final_genes,
        "intercept": logit_model_sm.params["const"],
        "gene_weights": gene_coefs["Coef"].to_dict(),
        # Youden-optimal cut-off on the test set — reused as the decision
        # threshold by predict_new_samples.py
        "youden_threshold": float(xgb_results["Cut-off"]),
    }
    joblib.dump(model_assets,
                os.path.join(timestamped_results_dir,
                             "Final_SCAP_Model_Assets.pkl"))
    print("Model assets saved to: Final_SCAP_Model_Assets.pkl")


if __name__ == "__main__":
    train_evaluate_final_panel(timestamped_results_dir)
