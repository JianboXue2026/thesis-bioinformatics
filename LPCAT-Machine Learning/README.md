# Clinical Data — DNN Classification

Deep neural network for binary classification on tabular clinical data.
Includes model training, performance evaluation, and comparative ROC analysis.

---

## Dependencies

| Package       | Version (tested) | Purpose                        |
| ------------- | ---------------- | ------------------------------ |
| Python        | >= 3.8           | Runtime                        |
| tensorflow    | >= 2.6           | Model building / training / TF ops |
| numpy         | —                | Array operations               |
| pandas        | —                | CSV I/O, data manipulation     |
| scikit-learn  | —                | Train-test split, ROC/AUC metrics |
| matplotlib    | —                | ROC curve plotting             |

### Install

```bash
pip install tensorflow numpy pandas scikit-learn matplotlib
```

---

## Project Structure

```
project/
│
├── data/                          # Example CSV datasets (full data archived at Peking University)
│   ├── MachineLearning.csv        # All 8 features present
│   ├── MachineLearning-nopct.csv  # PCT  → replaced with -1
│   ├── MachineLearning-nowbc.csv  # WBC  → replaced with -1
│   ├── MachineLearning-noesr.csv  # ESR  → replaced with -1
│   ├── MachineLearning-nocrp.csv  # CRP  → replaced with -1
│   ├── MachineLearning-noly.csv   # LY   → replaced with -1
│   ├── MachineLearning-nolpcat.csv# LPCAT1 → replaced with -1
│   └── MachineLearning-none.csv   # NE   → replaced with -1
│
├── models/                        # Saved Keras models (auto-generated)
│   └── <YYYYMMDD-HHMMSS>/         # One sub-folder per training run
│
├── results/roc/                   # ROC curve output images
│   └── <timestamp>_*.tiff
│
├── Machine_Learning.py            # Train a model from scratch
├── model_performance.py           # Evaluate a saved model + single ROC
├── ROCs.py                        # Compare ROC across dataset variants
├── Self_Defining_Function.py      # Shared prediction utilities
└── README.md
```

### CSV Format

All datasets share the same 27-column layout (zero-indexed):

| Col | Name        | Role     | Description                                |
|-----|-------------|----------|--------------------------------------------|
| 0   | `num`       | —        | Patient ID (not used as a feature)         |
| 1   | `SCAP`      | **Label**| Binary target: 0 = negative, 1 = positive  |
| 2   | `Death`     | —        | Mortality outcome                          |
| 3   | `Gender`    | —        | 0 = male, 1 = female                       |
| 4   | `Age`       | —        | Age in years                               |
| 5   | `BMI`       | —        | Body mass index                            |
| 6   | `Temperature`| —       | Body temperature (°C)                      |
| 7   | `RR`        | —        | Respiratory rate                           |
| 8   | `HR`        | —        | Heart rate                                 |
| 9   | `Conscious` | —        | 0 = alert, 1 = altered                     |
| 10  | `LPCAT1`    | **★**    | Lysophosphatidylcholine acyltransferase 1  |
| 11  | `WBC`       | **★**    | White blood cell count                     |
| 12  | `NE`        | **★**    | Neutrophil percentage                      |
| 13  | `LY`        | **★**    | Lymphocyte percentage                      |
| 14  | `NLR`       | **★**    | Neutrophil-to-lymphocyte ratio             |
| 15  | `Hb`        | —        | Hemoglobin                                 |
| 16  | `PLT`       | —        | Platelet count                             |
| 17  | `BUN`       | —        | Blood urea nitrogen                        |
| 18  | `Scr`       | —        | Serum creatinine                           |
| 19  | `ALB`       | —        | Albumin                                    |
| 20  | `CK`        | —        | Creatine kinase                            |
| 21  | `ESR`       | **★**    | Erythrocyte sedimentation rate             |
| 22  | `CRP`       | **★**    | C-reactive protein                         |
| 23  | `PCT`       | **★**    | Procalcitonin                              |
| 24  | `Glu`       | —        | Glucose                                    |
| 25  | `CURB65`    | —        | CURB-65 severity score (0–5)              |
| 26  | `PSI`       | —        | Pneumonia Severity Index                   |

**★** = feature columns used by the model (indices **10, 11, 12, 13, 14, 21, 22, 23**).

The label column (`SCAP`, index 1) contains `0` (negative) or `1` (positive).
Missing values are represented by `-1`.

#### Dataset Variants

Each dataset variant keeps the **same** 27-column structure but replaces one feature
column with `-1` across all rows. This allows testing the model's dependence on
each individual feature without changing the input shape.

| File                          | Dropped Feature | Column | Legend Label in ROCs.py |
|-------------------------------|-----------------|--------|-------------------------|
| `MachineLearning.csv`         | *(none)*        | —      | Prediction Model        |
| `MachineLearning-nopct.csv`   | PCT             | 23     | Drop-out PCT            |
| `MachineLearning-nowbc.csv`   | WBC             | 11     | Drop-out WBC            |
| `MachineLearning-noesr.csv`   | ESR             | 21     | Drop-out ESR            |
| `MachineLearning-nocrp.csv`   | CRP             | 22     | Drop-out CRP            |
| `MachineLearning-noly.csv`    | LY              | 13     | Drop-out LY%            |
| `MachineLearning-nolpcat.csv` | LPCAT1          | 10     | Drop-out LPCAT          |
| `MachineLearning-none.csv`    | NE              | 12     | Drop-out NE%            |

---

## Scripts

### 1. `Machine_Learning.py` — Train

Builds a 10-layer DNN and trains on the full dataset.

```
python Machine_Learning.py
```

- Reads `data/MachineLearning.csv`
- Splits into train / validation / test (60 / 20 / 20, stratified)
- Trains for 500 epochs with Adam + SparseCategoricalCrossentropy
- Saves model to `models/<timestamp>/`
- Prints test accuracy, loss, and confusion matrix

### 2. `model_performance.py` — Evaluate

Loads a saved model and evaluates it on a single dataset variant.

```
python model_performance.py
```

- Set `MODEL_TIMESTAMP` in the config to match a saved model folder
- Default dataset: `data/MachineLearning-nopct.csv`
- Outputs accuracy, loss, confusion matrix, AUC, and a ROC curve to `results/roc/`

### 3. `ROCs.py` — Compare ROC

Runs the same saved model against all 8 dataset variants and plots their ROC
curves on a single figure.

```
python ROCs.py
```

- Expects all 8 CSV files in `data/`
- Generates `results/roc/<timestamp>_ROCs.tiff`

### 4. `Self_Defining_Function.py` — Utilities

Shared helper functions used by the other scripts. Not meant to be run directly.

| Function           | Used by                  | Purpose                               |
|--------------------|--------------------------|---------------------------------------|
| `prediction()`     | `Machine_Learning.py`    | Softmax wrapper, returns class labels |
| `prediction_v30()` | `model_performance.py`   | Prediction + metrics + ROC plot       |
| `random_test()`    | (standalone utility)     | Random subset sampling                |
| `auto_log_to_csv()`| (standalone utility)     | Append experiment results to CSV log  |

---

## Workflow

1. **Train** — run `Machine_Learning.py`, note the model timestamp.
2. **Evaluate** — copy the timestamp into `model_performance.py` config, then run it.
3. **Compare** — copy the timestamp into `ROCs.py`, then run it to see all ROC curves.

---

## Notes

- The model uses `from_logits=True`; the final Dense layer has no activation.
  Prediction helpers (`Self_Defining_Function`) prepend a Softmax layer automatically.
- All random seeds are pinned (`RANDOM_STATE = 1`) for reproducibility.
- Edit the `FEATURE_COLUMNS` list in each script if your CSV layout differs.
