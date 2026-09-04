#!/usr/bin/env python3
"""Compare system and class balanced calibration on an external system."""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from calibration_utils import TRAIN_PATH, RANDOM_SEED, canonicalize


ROOT = Path(__file__).resolve().parents[1]
DATASET = sys.argv[1]
PREDICTORS = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
C_GRID = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]


def sample_weights(data, mode):
    if mode == "system":
        weights = data["dataset"].map(1 / data["dataset"].value_counts())
    elif mode == "system_class":
        counts = data.groupby(["dataset", "interaction_label"]).size()
        class_counts = data.groupby("dataset")["interaction_label"].nunique()
        weights = pd.Series(
            [
                1 / (class_counts[dataset] * counts[dataset, label])
                for dataset, label in zip(data["dataset"], data["interaction_label"])
            ],
            index=data.index,
        )
    else:
        raise ValueError(f"Unknown weighting mode: {mode}")
    return (weights / weights.mean()).to_numpy(float)


def make_model(c_value, fit_intercept=True):
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2",
            solver="liblinear",
            C=c_value,
            fit_intercept=fit_intercept,
            max_iter=5000,
            random_state=RANDOM_SEED,
        )),
    ])


def fit_model(data, mode, c_value, fit_intercept):
    model = make_model(c_value, fit_intercept)
    model.fit(
        data[PREDICTORS],
        data["interaction_label"],
        model__sample_weight=sample_weights(data, mode),
    )
    return model


def logo_predictions(data, mode, c_value, fit_intercept):
    predictions = pd.Series(index=data.index, dtype=float)
    for held_out in sorted(data["dataset"].unique()):
        training = data[data["dataset"] != held_out]
        model = fit_model(training, mode, c_value, fit_intercept)
        mask = data["dataset"] == held_out
        predictions.loc[mask] = model.predict_proba(data.loc[mask, PREDICTORS])[:, 1]
    return predictions


def macro_ap_lift(data, scores):
    lifts = []
    scored = data.assign(_score=scores)
    for _, group in scored.groupby("dataset"):
        if group["interaction_label"].nunique() == 2:
            prevalence = group["interaction_label"].mean()
            lifts.append(average_precision_score(group["interaction_label"], group["_score"]) / prevalence)
    return float(np.mean(lifts))


def select_c(data, mode, fit_intercept):
    candidates = []
    for c_value in C_GRID:
        scores = logo_predictions(data, mode, c_value, fit_intercept)
        candidates.append((macro_ap_lift(data, scores), -c_value, c_value))
    return max(candidates)[2], max(candidates)[0]


def binary_metrics(labels, selected):
    labels = np.asarray(labels, dtype=int)
    selected = np.asarray(selected, dtype=int)
    true_positive = int(((labels == 1) & (selected == 1)).sum())
    false_positive = int(((labels == 0) & (selected == 1)).sum())
    false_negative = int(((labels == 1) & (selected == 0)).sum())
    precision = true_positive / (true_positive + false_positive) if true_positive + false_positive else np.nan
    recall = true_positive / (true_positive + false_negative) if true_positive + false_negative else np.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else np.nan
    return precision, recall, f1


def evaluate(name, data, labels, scores, c_value, training_lift, model):
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    budget = int(labels.sum())
    selected = scores >= 0.5
    precision, recall, f1 = binary_metrics(labels, selected)
    order = np.argsort(-scores, kind="stable")
    density_selected = np.zeros(len(labels), dtype=int)
    density_selected[order[:budget]] = 1
    density_precision, density_recall, density_f1 = binary_metrics(labels, density_selected)
    fitted = model.named_steps["model"]
    return {
        "method": name,
        "tested_pairs": len(labels),
        "positive_pairs": budget,
        "neutral_pairs": int((labels == 0).sum()),
        "selected_C": c_value,
        "training_macro_ap_lift": training_lift,
        "average_precision": average_precision_score(labels, scores),
        "auroc": roc_auc_score(labels, scores),
        "threshold": 0.5,
        "selected_edges": int(selected.sum()),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "density_matched_edges": budget,
        "density_matched_true_positives": int(labels[density_selected == 1].sum()),
        "density_matched_precision": density_precision,
        "density_matched_recall": density_recall,
        "density_matched_f1": density_f1,
        "intercept": float(fitted.intercept_[0]) if fitted.fit_intercept else 0.0,
        **{
            f"coefficient_{predictor}": float(value)
            for predictor, value in zip(PREDICTORS, fitted.coef_[0])
        },
    }


def main():
    train = canonicalize(pd.read_csv(TRAIN_PATH))
    truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_tested_pairs.csv"))
    features = canonicalize(pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_pair_features.csv"))
    test = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one")
    labels = test["interaction_label"].astype(int).to_numpy()

    specifications = [
        ("system_weighted", "system", True, False),
        ("system_and_class_balanced", "system_class", True, False),
        ("system_and_class_balanced_no_intercept", "system_class", False, False),
        ("system_and_class_balanced_excluding_single_class", "system_class", True, True),
    ]
    results = []
    predictions = test[["dataset", "taxon_1", "taxon_2", "interaction_label", "tested_status"]].copy()
    for name, mode, fit_intercept, exclude_single_class in specifications:
        training = train.copy()
        if exclude_single_class:
            keep = training.groupby("dataset")["interaction_label"].transform("nunique") == 2
            training = training[keep].copy()
        c_value, training_lift = select_c(training, mode, fit_intercept)
        model = fit_model(training, mode, c_value, fit_intercept)
        scores = model.predict_proba(test[PREDICTORS])[:, 1]
        predictions[name] = scores
        results.append(evaluate(name, test, labels, scores, c_value, training_lift, model))

    result = pd.DataFrame(results)
    result.to_csv(ROOT / "results" / f"{DATASET}_balanced_calibration_performance.csv", index=False)
    predictions.to_csv(ROOT / "analysis_data" / f"{DATASET}_balanced_calibration_predictions.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
