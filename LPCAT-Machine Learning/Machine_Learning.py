"""
Machine Learning Model Training Script
--------------------------------------
Builds and trains a deep neural network on tabular clinical data.
The model is saved with a timestamp for later evaluation and prediction.
"""

import numpy as np
import pandas as pd
import tensorflow as tf
from datetime import datetime
from sklearn.model_selection import train_test_split
import Self_Defining_Function

# ============================================================
# Configuration
# ============================================================
CSV_DATA_PATH = "data/MachineLearning.csv"
MODEL_SAVE_DIR = "models/"

RANDOM_STATE = 1
TEST_SIZE = 0.2
VALID_SIZE = 0.2
EPOCHS = 500
BATCH_SIZE = 4

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
# Model Architecture
# ============================================================
model = tf.keras.Sequential([
    tf.keras.layers.Dense(8, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(128, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(128, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(64, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(64, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(32, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(32, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(16, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(8, activation='sigmoid', kernel_initializer='he_normal'),
    tf.keras.layers.Dense(2, kernel_regularizer=tf.keras.regularizers.L2(0.01),
                          kernel_initializer='he_normal'),
])

model.compile(
    optimizer=tf.keras.optimizers.Adam(),
    loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    metrics=['sparse_categorical_accuracy']
)

model.summary()

# ============================================================
# Timestamp for Model Saving
# ============================================================
timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

# ============================================================
# Training
# ============================================================
history = model.fit(
    X_train, y_train,
    epochs=EPOCHS,
    batch_size=BATCH_SIZE,
    validation_data=(X_valid, y_valid)
)

# ============================================================
# Save Model
# ============================================================
model_save_path = MODEL_SAVE_DIR + timestamp
model.save(model_save_path)
print(f"Model saved to: {model_save_path}")

# ============================================================
# Evaluation
# ============================================================
test_loss, test_acc = model.evaluate(X_test, y_test)
print(f"Test accuracy: {test_acc:.4f}")
print(f"Test loss:    {test_loss:.4f}")

# ============================================================
# Prediction & Confusion Matrix
# ============================================================
pred_labels = Self_Defining_Function.prediction(model, X_test)
print(f"Predicted labels: {pred_labels}")

cf_matrix = tf.math.confusion_matrix(y_test, pred_labels, num_classes=2)
print(f"Confusion Matrix:\n{cf_matrix}")
