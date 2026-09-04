#!/usr/bin/env python3
"""Evaluate compact supervised calibration models on the external Schäfer test."""

from pathlib import Path
import json
import math

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score, brier_score_loss

from calibration_utils import (
    TRAIN_PATH, TEST_DATASET, canonicalize, fit, logo_predictions,
    macro_f1, make_model, mean_ap_lift,
)


ROOT = Path(__file__).resolve().parents[1]
METHOD_SELECTION_COLUMNS = [
    "pln_selected", "glmnet_selected", "spieceasi_selected",
    "spring_selected", "sparcc_selected",
]
MODEL_SPECS = {
    "score_model": [
        "pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile",
    ],
    "observability_model": ["prevalence_min", "joint_prevalence"],
    "compact_combined_model": [
        "pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile",
        "method_agreement", "prevalence_min", "joint_prevalence",
    ],
}


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


def evaluate(data, score, selected, model_name, scope, threshold):
    y = data.interaction_label.astype(int).to_numpy()
    score = np.asarray(score, dtype=float)
    selected = np.asarray(selected, dtype=int)
    tp, fp, fn, precision, recall, f1 = binary_metrics(y, selected)
    budget = int(y.sum())
    order = np.lexsort((data.taxon_2, data.taxon_1, -score))
    density = np.zeros(len(data), dtype=int)
    density[order[:budget]] = 1
    dtp, dfp, dfn, dp, dr, df = binary_metrics(y, density)
    return {
        "scope": scope, "method": model_name, "tested_pairs": len(data),
        "positive_pairs": budget, "neutral_pairs": int((y == 0).sum()),
        "positive_prevalence": float(y.mean()),
        "average_precision": average_precision_score(y, score),
        "auroc": roc_auc_score(y, score), "brier_score": brier_score_loss(y, score),
        "mean_predicted_probability": float(score.mean()), "threshold": threshold,
        "native_selected_edges": int(selected.sum()), "true_positives": tp,
        "false_positives": fp, "false_negatives": fn,
        "native_precision": precision, "native_recall": recall, "native_f1": f1,
        "density_matched_edges": budget, "density_matched_true_positives": dtp,
        "density_matched_precision": dp, "density_matched_recall": dr,
        "density_matched_f1": df,
    }


def main():
    train = canonicalize(pd.read_csv(TRAIN_PATH))
    truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / f"{TEST_DATASET}_tested_pairs.csv"))
    features = canonicalize(pd.read_csv(ROOT / "analysis_data" / f"{TEST_DATASET}_pair_features.csv"))
    test = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one")
    for frame in (train, test):
        frame["method_agreement"] = frame[METHOD_SELECTION_COLUMNS].sum(axis=1)

    eligible = pd.read_csv(ROOT / "analysis_data" / f"{TEST_DATASET}_conet_tested_pair_scores.csv")
    eligible_pairs = set(zip(eligible.taxon_1.astype(str), eligible.taxon_2.astype(str)))
    conet_mask = test.apply(lambda row: (row.taxon_1, row.taxon_2) in eligible_pairs, axis=1)

    results = []
    tuning_rows = []
    coefficient_rows = []
    predictions = test[["dataset", "taxon_1", "taxon_2", "interaction_label", "tested_status"]].copy()
    metadata = {}

    for model_name, predictors in MODEL_SPECS.items():
        best = None
        for c_value in [0.001, 0.01, 0.1, 1.0, 10.0]:
            oof = logo_predictions(train, predictors, c_value, 0.0)
            objective = mean_ap_lift(train, oof)
            tuning_rows.append({
                "model": model_name, "C": c_value, "mean_ap_lift": objective,
            })
            candidate = (objective, -c_value, c_value, oof)
            if best is None or candidate[:2] > best[:2]:
                best = candidate

        _, _, c_value, oof = best
        thresholds = np.unique(np.r_[0.01, np.linspace(0.05, 0.95, 181), 0.99])
        threshold_values = [(macro_f1(train, oof, value), value) for value in thresholds]
        training_f1, threshold = max(threshold_values, key=lambda item: (item[0], item[1]))
        model = fit(make_model(c_value, 0.0), train[predictors], train.interaction_label, train.dataset)
        score = model.predict_proba(test[predictors])[:, 1]
        selected = (score >= threshold).astype(int)
        predictions[f"{model_name}_score"] = score
        predictions[f"{model_name}_selected"] = selected

        results.append(evaluate(test, score, selected, model_name, "all_schafer_tested_pairs", threshold))
        results.append(evaluate(
            test.loc[conet_mask], score[conet_mask.to_numpy()], selected[conet_mask.to_numpy()],
            model_name, "conet_eligible_pairs", threshold,
        ))
        for predictor, coefficient in zip(predictors, model.named_steps["model"].coef_[0]):
            coefficient_rows.append({
                "model": model_name, "predictor": predictor,
                "standardized_coefficient": coefficient,
            })
        metadata[model_name] = {
            "predictors": predictors, "selected_C": c_value,
            "selected_threshold": threshold, "training_macro_f1": training_f1,
            "test_labels_used_for_fitting_or_tuning": False,
        }

    rank_baselines = {
        "mean_three_score_rank": test[
            ["pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile"]
        ].mean(axis=1).to_numpy(),
        "mean_pln_sparcc_rank": test[
            ["pln_score_percentile", "sparcc_score_percentile"]
        ].mean(axis=1).to_numpy(),
    }
    for baseline_name, score in rank_baselines.items():
        no_native_graph = np.zeros(len(test), dtype=int)
        results.append(evaluate(
            test, score, no_native_graph, baseline_name,
            "all_schafer_tested_pairs", math.nan,
        ))
        results.append(evaluate(
            test.loc[conet_mask], score[conet_mask.to_numpy()],
            no_native_graph[conet_mask.to_numpy()], baseline_name,
            "conet_eligible_pairs", math.nan,
        ))

    result_frame = pd.DataFrame(results)
    result_frame.to_csv(ROOT / "results" / "compact_calibration_experiments_schafer.csv", index=False)
    pd.DataFrame(tuning_rows).to_csv(ROOT / "results" / "compact_calibration_tuning.csv", index=False)
    pd.DataFrame(coefficient_rows).to_csv(ROOT / "results" / "compact_calibration_coefficients.csv", index=False)
    predictions.to_csv(ROOT / "analysis_data" / f"{TEST_DATASET}_compact_calibration_predictions.csv", index=False)
    with open(ROOT / "results" / "compact_calibration_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)

    print(result_frame.to_string(index=False))
    print("\nCoefficients")
    print(pd.DataFrame(coefficient_rows).to_string(index=False))


if __name__ == "__main__":
    main()
