"""
Model Performance Evaluation Script
------------------------------------
Loads a previously trained model, evaluates it on test data,
and generates predictions, confusion matrix, and ROC curve.
"""

import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
import Self_Defining_Function

# ============================================================
# Configuration
# ============================================================
CSV_DATA_PATH = "data/MachineLearning-nopct.csv"
MODEL_DIR = "models/"
ROC_OUTPUT_DIR = "results/roc/"
# Set the model timestamp to load (generated during training)
MODEL_TIMESTAMP = "20250902-005017"

RANDOM_STATE = 1
TEST_SIZE = 0.2
VALID_SIZE = 0.2

# Feature column indices (0-based) in the CSV
FEATURE_COLUMNS = [10, 11, 12, 13, 14, 21, 22, 23]
# Label column index (0-based)
LABEL_COLUMN = [1]

# ============================================================
# Data Loading
# ============================================================
df = pd.read_csv(CSV_DATA_PATH, sep=',')
print(f"Dataset loaded — shape: {df.shape}")

X = df.iloc[:, FEATURE_COLUMNS]
y = df.iloc[:, LABEL_COLUMN]

# ============================================================
# Train / Test Split
# ============================================================
X_train_full, X_test, y_train_full, y_test = train_test_split(
    np.array(X), np.array(y),
    test_size=TEST_SIZE,
    random_state=RANDOM_STATE,
    stratify=y
)

# ============================================================
# Train / Validation Split
# ============================================================
X_train, X_valid, y_train, y_valid = train_test_split(
    X_train_full, y_train_full,
    test_size=VALID_SIZE,
    random_state=RANDOM_STATE,
    stratify=y_train_full
)

print(f"Train: {X_train.shape}  |  Validation: {X_valid.shape}  |  Test: {X_test.shape}")

# ============================================================
# Load Saved Model
# ============================================================
model_path = MODEL_DIR + MODEL_TIMESTAMP
model = tf.keras.models.load_model(model_path)
print(f"Model loaded from: {model_path}")

# ============================================================
# Evaluation
# ============================================================
test_loss, test_acc = model.evaluate(X_test, y_test)
print(f"Test accuracy: {test_acc:.4f}")
print(f"Test loss:    {test_loss:.4f}")

# ============================================================
# Prediction, Confusion Matrix & ROC Curve
# ============================================================
roc_path = ROC_OUTPUT_DIR + MODEL_TIMESTAMP + "_nopct.tiff"

pred_labels, auc, cf_matrix = Self_Defining_Function.prediction_v30(
    origin_model=model,
    test_data=X_test,
    real_label=y_test,
    roc_pic_path=roc_path
)

print(f"Predicted labels: {pred_labels}")
print(f"Confusion Matrix:\n{cf_matrix}")
print(f"AUC: {auc:.4f}")
