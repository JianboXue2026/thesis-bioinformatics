"""
ROC Curve Generation Script
----------------------------
Evaluates a saved model across multiple datasets (different feature subsets)
and plots comparative ROC curves for each variant.
"""

import numpy as np
import pandas as pd
import tensorflow as tf
import sklearn.metrics
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

# ============================================================
# Configuration
# ============================================================
# List of dataset variants to evaluate
DATASET_NAMES = [
    "data/MachineLearning.csv",
    "data/MachineLearning-nopct.csv",
    "data/MachineLearning-nowbc.csv",
    "data/MachineLearning-noesr.csv",
    "data/MachineLearning-nocrp.csv",
    "data/MachineLearning-noly.csv",
    "data/MachineLearning-nolpcat.csv",
    "data/MachineLearning-none.csv",
]

MODEL_DIR = "models/"
ROC_OUTPUT_DIR = "results/roc/"
# Model timestamp to load
MODEL_TIMESTAMP = "20250902-005017"

RANDOM_STATE = 1
TEST_SIZE = 0.2
VALID_SIZE = 0.2

# Feature column indices (0-based)
FEATURE_COLUMNS = [10, 11, 12, 13, 14, 21, 22, 23]
LABEL_COLUMN = [1]

# ============================================================
# Load Model
# ============================================================
model_path = MODEL_DIR + MODEL_TIMESTAMP
model = tf.keras.models.load_model(model_path)
print(f"Model loaded from: {model_path}")

# ============================================================
# Evaluate on Each Dataset Variant
# ============================================================
roc_data_list = []  # Stores [fpr, tpr, auc] for each dataset

for dataset_name in DATASET_NAMES:
    print(f"\n{'='*60}")
    print(f"Processing: {dataset_name}")

    # --- Load and split data ---
    df = pd.read_csv(dataset_name, sep=',')
    X = df.iloc[:, FEATURE_COLUMNS]
    y = df.iloc[:, LABEL_COLUMN]

    X_train_full, X_test, y_train_full, y_test = train_test_split(
        np.array(X), np.array(y),
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y
    )

    # --- Evaluate model ---
    test_loss, test_acc = model.evaluate(X_test, y_test, verbose=0)
    print(f"Test accuracy: {test_acc:.4f}  |  Test loss: {test_loss:.4f}")

    # --- Generate predictions ---
    pred = model.predict(X_test, verbose=0)

    pred_label_list = [np.argmax(item) for item in pred]

    # --- Confusion Matrix ---
    cf_matrix = tf.math.confusion_matrix(y_test, pred_label_list, num_classes=2)
    print(f"Confusion Matrix:\n{cf_matrix}")

    # --- ROC Curve Data ---
    pred_prob_list = [item[1] for item in pred]  # Probability of positive class (index 1)

    fpr, tpr, _ = sklearn.metrics.roc_curve(y_test, pred_prob_list, pos_label=1)
    auc = sklearn.metrics.auc(fpr, tpr)
    print(f"AUC: {auc:.4f}")

    roc_data_list.append([fpr, tpr, auc])

# ============================================================
# Plot Combined ROC Curves
# ============================================================
# Legend labels for each dataset variant
legend_names = [
    "Prediction Model",
    "Drop-out PCT",
    "Drop-out WBC",
    "Drop-out ESR",
    "Drop-out CRP",
    "Drop-out LY%",
    "Drop-out LPCAT",
    "Drop-out NE%",
]

# Color scheme for each curve
colors = ['red', 'blue', 'orange', 'green', 'purple', 'brown', 'olive', 'cyan']

fig = plt.figure()
lw = 2

for idx, (fpr, tpr, auc) in enumerate(roc_data_list):
    plt.plot(fpr, tpr, color=colors[idx], lw=lw,
             label=f"{legend_names[idx]} (AUC = {auc:.2f})")

# Reference line (random classifier)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('1 - Specificity (False Positive Rate)')
plt.ylabel('Sensitivity (True Positive Rate)')
plt.legend(loc="lower right")

output_path = ROC_OUTPUT_DIR + MODEL_TIMESTAMP + "_ROCs.tiff"
fig.savefig(output_path, dpi=600, format='tiff')
print(f"\nROC curves saved to: {output_path}")
plt.show()
