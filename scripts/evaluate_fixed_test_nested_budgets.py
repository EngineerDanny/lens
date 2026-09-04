#!/usr/bin/env python3
"""Evaluate nested calibration budgets on a fixed test set in each repeat."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from evaluate_direct_pair_features import FEATURE_SETS, load_analysis, score_model
from evaluate_within_system_calibration import ROOT, SEED, UNSUPERVISED


SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
BUDGET_FRACTIONS = [0.20, 0.40, 0.60, 0.80]
TEST_FRACTION = 0.20
REPEATS = 200
METHODS = {
    **UNSUPERVISED,
    "SPIEC-EASI": "spieceasi_score_percentile",
}


def draw_fixed_test(frame: pd.DataFrame, pair_ids: np.ndarray,
                    test_size: int, rng: np.random.Generator) -> set[str]:
    """Draw an evaluation set with both outcomes; labels only validate the split."""
    for _ in range(10_000):
        chosen = set(rng.choice(pair_ids, size=test_size, replace=False))
        test = frame[frame.pair_id.isin(chosen)]
        if test.interaction_label.nunique() == 2:
            return chosen
    raise RuntimeError("Could not draw a test set containing both outcomes")


def run_system(system: str) -> tuple[list[dict], list[dict]]:
    frame = load_analysis(system)
    pair_ids = frame.pair_id.drop_duplicates().to_numpy()
    n_pairs = len(pair_ids)
    test_size = max(2, int(np.ceil(TEST_FRACTION * n_pairs)))
    budget_sizes = {
        fraction: int(np.ceil(fraction * n_pairs)) for fraction in BUDGET_FRACTIONS
    }
    rows, failures = [], []

    for repeat in range(REPEATS):
        rng = np.random.default_rng(SEED + 1_000_000 * SYSTEMS.index(system) + repeat)
        try:
            test_ids = draw_fixed_test(frame, pair_ids, test_size, rng)
        except RuntimeError as error:
            failures.append({"analysis_set": system, "repeat": repeat, "reason": str(error)})
            continue
        pool = np.array([pair for pair in pair_ids if pair not in test_ids], dtype=object)
        ordered_pool = rng.permutation(pool)
        test = frame[frame.pair_id.isin(test_ids)]
        y_test = test.interaction_label.to_numpy()

        # These scores and test pairs are identical at every budget in this repeat.
        for method, column in METHODS.items():
            score = test[column].to_numpy()
            rows.append({
                "analysis_set": system,
                "repeat": repeat,
                "budget_fraction": 0.0,
                "budget_pairs": 0,
                "test_pairs": len(test_ids),
                "test_positives": int(y_test.sum()),
                "method": method,
                "auprc": average_precision_score(y_test, score),
                "auroc": roc_auc_score(y_test, score),
                "model_fitted": False,
            })

        for fraction in BUDGET_FRACTIONS:
            budget = min(budget_sizes[fraction], len(ordered_pool))
            calibration_ids = set(ordered_pool[:budget])
            calibration = frame[frame.pair_id.isin(calibration_ids)]
            _, score, fitted = score_model(calibration, test, FEATURE_SETS["combined"])
            rows.append({
                "analysis_set": system,
                "repeat": repeat,
                "budget_fraction": fraction,
                "budget_pairs": budget,
                "test_pairs": len(test_ids),
                "test_positives": int(y_test.sum()),
                "method": "Supervised combined",
                "auprc": average_precision_score(y_test, score),
                "auroc": roc_auc_score(y_test, score),
                "model_fitted": fitted,
            })
    return rows, failures


def main() -> None:
    rows, failures = [], []
    for system in SYSTEMS:
        system_rows, system_failures = run_system(system)
        rows.extend(system_rows)
        failures.extend(system_failures)
    results = pd.DataFrame(rows)
    summary = results.groupby(
        ["analysis_set", "method", "budget_fraction"], as_index=False
    ).agg(
        valid_repeats=("repeat", "size"),
        budget_pairs=("budget_pairs", "mean"),
        test_pairs=("test_pairs", "mean"),
        test_positives=("test_positives", "mean"),
        model_fit_rate=("model_fitted", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
    )
    results.to_csv(ROOT / "results" / "fixed_test_nested_budget_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "fixed_test_nested_budget_summary.csv", index=False)
    pd.DataFrame(failures, columns=["analysis_set", "repeat", "reason"]).to_csv(
        ROOT / "results" / "fixed_test_nested_budget_failures.csv", index=False
    )
    with open(ROOT / "results" / "fixed_test_nested_budget_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "repeats": REPEATS,
            "test_fraction": TEST_FRACTION,
            "budget_fractions_of_all_tested_pair_identifiers": BUDGET_FRACTIONS,
            "test_set_fixed_across_budgets_within_repeat": True,
            "calibration_sets_nested_within_repeat": True,
            "calibration_pair_sampling_uses_outcomes": False,
            "test_split_outcomes_used_only_to_require_both_classes": True,
            "model": "fixed L2 logistic regression, C=0.01, balanced class weights",
            "combined_features": FEATURE_SETS["combined"],
            "unsupervised_methods": list(METHODS),
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
