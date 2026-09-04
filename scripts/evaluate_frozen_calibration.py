#!/usr/bin/env python3
"""Apply the pre-specified three-score calibration model to a new system."""

from pathlib import Path
import math
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from calibration_utils import TRAIN_PATH, canonicalize, fit, make_model


ROOT = Path(__file__).resolve().parents[1]
DATASET = sys.argv[1]
PREDICTORS = ["pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile"]
FROZEN_C = 1.0
FROZEN_THRESHOLD = 0.66


def binary_metrics(y, selected):
    y = np.asarray(y, dtype=int)
    selected = np.asarray(selected, dtype=int)
    tp = int(((y == 1) & (selected == 1)).sum())
    fp = int(((y == 0) & (selected == 1)).sum())
    fn = int(((y == 1) & (selected == 0)).sum())
    precision = tp / (tp + fp) if tp + fp else math.nan
    recall = tp / (tp + fn) if tp + fn else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else math.nan
    return tp, fp, fn, precision, recall, f1


def evaluate(name, y, score, selected):
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    selected = np.asarray(selected, dtype=int)
    tp, fp, fn, precision, recall, f1 = binary_metrics(y, selected)
    budget = int(y.sum())
    order = np.argsort(-score, kind="stable")
    density = np.zeros(len(y), dtype=int)
    density[order[:budget]] = 1
    dtp, _, _, dp, dr, df = binary_metrics(y, density)
    return {
        "method": name,
        "tested_pairs": len(y),
        "positive_pairs": budget,
        "neutral_pairs": int((y == 0).sum()),
        "average_precision": average_precision_score(y, score),
        "auroc": roc_auc_score(y, score),
        "brier_score": brier_score_loss(y, np.clip(score, 0, 1)),
        "selected_edges": int(selected.sum()),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "density_matched_edges": budget,
        "density_matched_true_positives": dtp,
        "density_matched_precision": dp,
        "density_matched_recall": dr,
        "density_matched_f1": df,
    }


def main():
    train = canonicalize(pd.read_csv(TRAIN_PATH))
    truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_tested_pairs.csv"))
    features = canonicalize(pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_pair_features.csv"))
    test = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one")
    model = fit(make_model(FROZEN_C), train[PREDICTORS], train.interaction_label, train.dataset)
    calibrated = model.predict_proba(test[PREDICTORS])[:, 1]
    mean_rank = test[PREDICTORS].mean(axis=1).to_numpy()
    y = test.interaction_label.astype(int).to_numpy()
    rows = [
        evaluate("frozen_three_score_calibration", y, calibrated, calibrated >= FROZEN_THRESHOLD),
        evaluate("mean_three_score_rank", y, mean_rank, np.zeros(len(y), dtype=int)),
    ]
    predictions = test[["dataset", "taxon_1", "taxon_2", "interaction_label", "tested_status"]].copy()
    predictions["calibrated_probability"] = calibrated
    predictions["frozen_selected"] = (calibrated >= FROZEN_THRESHOLD).astype(int)
    predictions["mean_three_score_rank"] = mean_rank
    predictions.to_csv(ROOT / "analysis_data" / f"{DATASET}_frozen_calibration_predictions.csv", index=False)
    result = pd.DataFrame(rows)
    result.to_csv(ROOT / "results" / f"{DATASET}_frozen_calibration_performance.csv", index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
