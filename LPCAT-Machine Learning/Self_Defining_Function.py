"""
Custom Utility Functions
------------------------
Project-specific helper functions for model prediction, evaluation,
random sampling, and experiment logging.
"""

import numpy as np
import pandas as pd
import tensorflow as tf
import sklearn.metrics
import matplotlib.pyplot as plt
import random

# ============================================================
# Model Prediction
# ============================================================


def prediction(origin_model, test_data):
    """
    Generate class predictions from a model that outputs logits.

    A Softmax layer is prepended to convert raw logits into
    class probabilities, then argmax selects the predicted class.

    Args:
        origin_model: A trained tf.keras model (logit output).
        test_data:   Input features as a numpy array.

    Returns:
        List of predicted class labels (integers).
    """
    model = tf.keras.Sequential()
    model.add(origin_model)
    model.add(tf.keras.layers.Softmax())

    pred = model.predict(test_data, verbose=0)
    pred_labels = [np.argmax(item) for item in pred]

    return pred_labels


def prediction_v30(origin_model_path=None, origin_model=None,
                   test_data=None, real_label=None, roc_pic_path=None):
    """
    Extended prediction with evaluation metrics and ROC curve generation.

    Supports loading a model either from a file path or by passing
    an already-loaded model object directly.

    Computes:
        - Predicted class labels
        - Confusion matrix
        - ROC curve and AUC
        - Saves the ROC plot to disk

    Args:
        origin_model_path: Path to a saved tf.keras model (optional).
        origin_model:      An already-loaded tf.keras model (optional).
        test_data:         Input features for prediction.
        real_label:        Ground-truth labels.
        roc_pic_path:      File path to save the ROC curve image.

    Returns:
        Tuple of (predicted_labels, auc, confusion_matrix).
    """
    # --- Load or use provided model ---
    if origin_model is None:
        model = tf.keras.models.load_model(origin_model_path)
    elif origin_model_path is None:
        model = origin_model
    else:
        raise ValueError("Provide either origin_model_path or origin_model, not both.")

    # --- Add Softmax layer for probability output ---
    pred_model = tf.keras.Sequential()
    pred_model.add(model)
    pred_model.add(tf.keras.layers.Softmax())

    # --- Predict ---
    pred = pred_model.predict(test_data, verbose=0)
    pred_labels = [np.argmax(item) for item in pred]
    print(f"Predictions: {pred_labels}")

    # --- Confusion Matrix ---
    cf_matrix = tf.math.confusion_matrix(real_label, pred_labels, num_classes=2)
    print(f"Confusion Matrix:\n{cf_matrix}")

    # --- ROC Curve & AUC ---
    pred_probs = [item[1] for item in pred]  # Probability for positive class (index 1)

    fpr, tpr, _ = sklearn.metrics.roc_curve(real_label, pred_probs, pos_label=1)
    auc = sklearn.metrics.auc(fpr, tpr)
    print(f"AUC = {auc:.4f}")

    # --- Plot & Save ROC Curve ---
    fig = plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='red', lw=lw,
             label=f'Model (AUC = {auc:.2f})')
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('1 - Specificity (False Positive Rate)')
    plt.ylabel('Sensitivity (True Positive Rate)')
    plt.legend(loc="lower right")
    fig.savefig(roc_pic_path, dpi=600, format='tiff')
    plt.show()

    return pred_labels, auc, cf_matrix


# ============================================================
# Utility Functions
# ============================================================


def random_test(random_num, test_img_all, test_label_all):
    """
    Randomly sample a subset of test data.

    Args:
        random_num:      Number of samples to draw.
        test_img_all:    Full set of test images/features.
        test_label_all:  Corresponding labels.

    Returns:
        Tuple of (sampled_features, sampled_labels) as numpy arrays.
    """
    sampled_features = []
    sampled_labels = []
    total = len(test_img_all)

    for _ in range(random_num):
        idx = random.randint(0, total - 1)
        sampled_features.append(test_img_all[idx])
        sampled_labels.append(test_label_all[idx])

    return np.array(sampled_features), sampled_labels


def auto_log_to_csv(model_name, csv_file, data_row):
    """
    Append a new row of experiment results to a CSV log file.

    Args:
        model_name: Identifier for the model run.
        csv_file:   Path to the CSV log file.
        data_row:   List of values to append, matching the column order.

    The CSV is expected to have the following columns:
        Logs_name, AHI, Position, Training_State, Epochs, Batch_Size,
        Training_Acc, Valid_Acc, Test_Acc, TN, FP, FN, TP, AUC,
        Dataset_Creation_Time, Model_Train_Time, Forward-Propagation_Time,
        Time_All
    """
    column_names = [
        'Logs_name', 'AHI', 'Position', 'Training_State', 'Epochs',
        'Batch_Size', 'Training_Acc', 'Valid_Acc', 'Test_Acc',
        'TN', 'FP', 'FN', 'TP', 'AUC',
        'Dataset_Creation_Time', 'Model_Train_Time',
        'Forward-Propagation_Time', 'Time_All'
    ]

    csv = pd.read_csv(csv_file, sep=',', index_col=None)
    new_row = pd.DataFrame([data_row], columns=column_names)
    updated = pd.concat([csv, new_row], axis=0, ignore_index=True)
    updated.to_csv(csv_file, index=None)
