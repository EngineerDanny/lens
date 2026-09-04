#!/usr/bin/env python3
"""Evaluate AUPRC-tuned sparse logistic regression in fixed outer folds."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from evaluate_pair_five_fold_cv import (
    BUDGET_FRACTIONS,
    REFERENCE_METHODS,
    SYSTEMS,
    load_system,
)
from evaluate_within_system_calibration import ROOT, SEED
from optimize_supervised_model import Candidate, fit_score, grouped_folds


CANDIDATES = [
    Candidate("centered_ridge", "lambda=0.01"),
    Candidate("centered_ridge", "lambda=0.1"),
    Candidate("centered_ridge", "lambda=1"),
    Candidate("centered_ridge", "lambda=10"),
]


def select_regularization(training: pd.DataFrame, seed: int):
    folds = grouped_folds(training, seed)
    if not folds:
        return CANDIDATES[0], []
    rows = []
    for candidate_index, candidate in enumerate(CANDIDATES):
        scores = []
        for fold_index, (inner_train, inner_valid) in enumerate(folds):
            y = inner_valid.interaction_label.to_numpy(int)
            if np.unique(y).size < 2:
                continue
            predicted = fit_score(
                candidate, inner_train, inner_valid,
                seed + 1000 * candidate_index + fold_index,
            )
            scores.append(average_precision_score(y, predicted))
        rows.append({
            "setting": candidate.setting,
            "inner_folds": len(scores),
            "mean_inner_auprc": float(np.mean(scores)) if scores else np.nan,
        })
    valid = [row for row in rows if np.isfinite(row["mean_inner_auprc"])]
    if not valid:
        return CANDIDATES[0], rows
    # Prefer stronger regularization when inner AUPRC ties.
    best = max(valid, key=lambda row: (row["mean_inner_auprc"],
                                       float(row["setting"].split("=")[1])))
    return Candidate("centered_ridge", best["setting"]), rows


def run_system(system: str):
    frame = load_system(system)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    _, onenet_column = REFERENCE_METHODS[system]
    output, tuning_output = [], []
    for fold, (pool_index, test_index) in enumerate(
        splitter.split(groups.pair_id, groups.group_label), start=1
    ):
        pool_ids = groups.iloc[pool_index].pair_id.to_numpy()
        test_ids = set(groups.iloc[test_index].pair_id)
        test = frame[frame.pair_id.isin(test_ids)]
        y_test = test.interaction_label.to_numpy(int)
        rng = np.random.default_rng(SEED + 100_000 * SYSTEMS.index(system) + fold)
        ordered_pool = rng.permutation(pool_ids)
        for fraction in BUDGET_FRACTIONS:
            budget = min(int(np.ceil(fraction * len(groups))), len(ordered_pool))
            train_ids = set(ordered_pool[:budget])
            training = frame[frame.pair_id.isin(train_ids)]
            seed = (SEED + 1_000_000 * SYSTEMS.index(system) +
                    10_000 * fold + int(100 * fraction))
            candidate, tuning = select_regularization(training, seed)
            predicted = fit_score(candidate, training, test, seed)
            output.append({
                "analysis_set": system,
                "fold": fold,
                "budget_fraction": fraction,
                "budget_pairs": budget,
                "test_pairs": len(test_ids),
                "test_positives": int(y_test.sum()),
                "selected_setting": candidate.setting,
                "auprc": average_precision_score(y_test, predicted),
                "auroc": roc_auc_score(y_test, predicted),
                "onenet_auprc": average_precision_score(y_test, test[onenet_column]),
            })
            for row in tuning:
                tuning_output.append({
                    "analysis_set": system, "fold": fold,
                    "budget_fraction": fraction, **row,
                })
    return output, tuning_output


def main():
    rows, tuning = [], []
    for system in SYSTEMS:
        print(f"Evaluating {system}", flush=True)
        system_rows, system_tuning = run_system(system)
        rows.extend(system_rows)
        tuning.extend(system_tuning)
    results = pd.DataFrame(rows)
    current = pd.read_csv(ROOT / "results" / "pair_five_fold_cv.csv")
    current = current[current.method.eq("Supervised combined")][
        ["analysis_set", "fold", "budget_fraction", "auprc"]
    ].rename(columns={"auprc": "current_auprc"})
    results = results.merge(
        current, on=["analysis_set", "fold", "budget_fraction"],
        validate="one_to_one",
    )
    summary = results.groupby(
        ["analysis_set", "budget_fraction"], as_index=False
    ).agg(
        budget_pairs=("budget_pairs", "max"), folds=("fold", "size"),
        sparse_auprc=("auprc", "mean"), sparse_sd=("auprc", "std"),
        current_auprc=("current_auprc", "mean"),
        onenet_auprc=("onenet_auprc", "mean"),
        sparse_auroc=("auroc", "mean"),
    )
    summary["improvement"] = summary.sparse_auprc - summary.current_auprc
    results.to_csv(ROOT / "results" / "optimized_sparse_logistic.csv", index=False)
    summary.to_csv(
        ROOT / "results" / "optimized_sparse_logistic_summary.csv", index=False
    )
    pd.DataFrame(tuning).to_csv(
        ROOT / "results" / "optimized_sparse_logistic_tuning.csv", index=False
    )
    with open(ROOT / "results" / "optimized_sparse_logistic_metadata.json", "w") as handle:
        json.dump({
            "model_family": "centered ridge logistic regression",
            "regularization_candidates_lambda": [0.01, 0.1, 1.0, 10.0],
            "coefficient_reference": "equal positive weights on PLNNetwork, Poisson GLMNet, and SparCC percentiles; zero weights on abundance features",
            "fallback": "fixed mean of the three score percentiles when either outcome class has fewer than three training pairs",
            "selection_metric": "mean inner-fold AUPRC",
            "outer_folds": 5, "inner_folds_maximum": 3,
            "test_labels_used_for_selection": False,
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
