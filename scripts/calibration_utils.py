"""Shared utilities for compact supervised calibration experiments."""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
TRAIN_PATH = ROOT / "training_data" / "all_labeled_pairs.csv"
TEST_DATASET = "schafer_phyllosphere_2022"
RANDOM_SEED = 20260820


def canonicalize(frame):
    ordered = frame[["taxon_1", "taxon_2"]].apply(
        lambda row: sorted(row.astype(str)), axis=1, result_type="expand"
    )
    frame = frame.copy()
    frame[["taxon_1", "taxon_2"]] = ordered.to_numpy()
    return frame


def study_weights(groups):
    counts = pd.Series(groups).value_counts()
    weights = pd.Series(groups).map(1 / counts).to_numpy(float)
    return weights / weights.mean()


def make_model(c_value, _unused_l1_ratio=0.0):
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2", solver="liblinear", C=c_value,
            max_iter=5000, random_state=RANDOM_SEED,
        )),
    ])


def fit(model, x, y, groups):
    model.fit(x, y, model__sample_weight=study_weights(groups))
    return model


def logo_predictions(train, predictors, c_value, l1_ratio=0.0):
    predictions = pd.Series(index=train.index, dtype=float)
    for held_out in sorted(train.dataset.unique()):
        fit_rows = train.dataset != held_out
        model = fit(
            make_model(c_value, l1_ratio), train.loc[fit_rows, predictors],
            train.loc[fit_rows, "interaction_label"], train.loc[fit_rows, "dataset"],
        )
        predictions.loc[~fit_rows] = model.predict_proba(train.loc[~fit_rows, predictors])[:, 1]
    return predictions


def mean_ap_lift(train, scores):
    values = []
    for _, group in train.assign(_score=scores).groupby("dataset"):
        y = group.interaction_label.astype(int)
        if y.nunique() == 2:
            values.append(average_precision_score(y, group._score) / y.mean())
    return float(np.mean(values))


def macro_f1(train, scores, threshold):
    values = []
    for _, group in train.assign(_score=scores).groupby("dataset"):
        y = group.interaction_label.astype(int)
        if y.nunique() < 2:
            continue
        selected = (group._score >= threshold).astype(int)
        tp = int(((y == 1) & (selected == 1)).sum())
        fp = int(((y == 0) & (selected == 1)).sum())
        fn = int(((y == 1) & (selected == 0)).sum())
        precision = tp / (tp + fp) if tp + fp else 0
        recall = tp / (tp + fn) if tp + fn else 0
        values.append(2 * precision * recall / (precision + recall) if precision + recall else 0)
    return float(np.mean(values))
