#!/usr/bin/env python3
"""Evaluate nested label budgets with a fixed held-out taxon panel."""

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
TEST_TAXON_FRACTION = 0.40
REPEATS = 200
METHODS = {
    **UNSUPERVISED,
    "SPIEC-EASI": "spieceasi_score_percentile",
}


def choose_partition(frame: pd.DataFrame, taxa: np.ndarray,
                     test_taxa_size: int, rng: np.random.Generator):
    """Choose disjoint test and calibration taxa with evaluable pair sets."""
    for _ in range(10_000):
        test_taxa = set(rng.choice(taxa, size=test_taxa_size, replace=False))
        calibration_taxa = set(taxa) - test_taxa
        test_ids = set(frame.loc[
            frame.taxon_1.isin(test_taxa) & frame.taxon_2.isin(test_taxa), "pair_id"
        ])
        pool_ids = set(frame.loc[
            frame.taxon_1.isin(calibration_taxa)
            & frame.taxon_2.isin(calibration_taxa), "pair_id"
        ])
        test = frame[frame.pair_id.isin(test_ids)]
        if test_ids and pool_ids and test.interaction_label.nunique() == 2:
            return test_taxa, calibration_taxa, test_ids, pool_ids
    raise RuntimeError("Could not form evaluable disjoint taxon panels")


def run_system(system: str) -> tuple[list[dict], list[dict]]:
    frame = load_analysis(system)
    taxa = np.array(sorted(set(frame.taxon_1) | set(frame.taxon_2)), dtype=object)
    test_taxa_size = min(
        max(3, int(round(TEST_TAXON_FRACTION * len(taxa)))), len(taxa) - 3
    )
    rows, failures = [], []

    for repeat in range(REPEATS):
        rng = np.random.default_rng(SEED + 2_000_000 * SYSTEMS.index(system) + repeat)
        try:
            test_taxa, calibration_taxa, test_ids, pool_ids = choose_partition(
                frame, taxa, test_taxa_size, rng
            )
        except RuntimeError as error:
            failures.append({"analysis_set": system, "repeat": repeat, "reason": str(error)})
            continue

        ordered_pool = rng.permutation(np.array(sorted(pool_ids), dtype=object))
        test = frame[frame.pair_id.isin(test_ids)]
        y_test = test.interaction_label.to_numpy()

        for method, column in METHODS.items():
            score = test[column].to_numpy()
            rows.append({
                "analysis_set": system,
                "repeat": repeat,
                "budget_fraction": 0.0,
                "budget_pairs": 0,
                "calibration_pool_pairs": len(pool_ids),
                "calibration_taxa": len(calibration_taxa),
                "test_taxa": len(test_taxa),
                "test_pairs": len(test_ids),
                "test_positives": int(y_test.sum()),
                "method": method,
                "auprc": average_precision_score(y_test, score),
                "auroc": roc_auc_score(y_test, score),
                "model_fitted": False,
            })

        for fraction in BUDGET_FRACTIONS:
            budget = min(max(1, int(np.ceil(fraction * len(pool_ids)))), len(pool_ids))
            calibration_ids = set(ordered_pool[:budget])
            calibration = frame[frame.pair_id.isin(calibration_ids)]
            _, score, fitted = score_model(calibration, test, FEATURE_SETS["combined"])
            rows.append({
                "analysis_set": system,
                "repeat": repeat,
                "budget_fraction": fraction,
                "budget_pairs": budget,
                "calibration_pool_pairs": len(pool_ids),
                "calibration_taxa": len(calibration_taxa),
                "test_taxa": len(test_taxa),
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
        mean_budget_pairs=("budget_pairs", "mean"),
        mean_calibration_pool_pairs=("calibration_pool_pairs", "mean"),
        calibration_taxa=("calibration_taxa", "mean"),
        test_taxa=("test_taxa", "mean"),
        mean_test_pairs=("test_pairs", "mean"),
        mean_test_positives=("test_positives", "mean"),
        model_fit_rate=("model_fitted", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
    )
    results.to_csv(ROOT / "results" / "fixed_taxon_panel_nested_budget_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "fixed_taxon_panel_nested_budget_summary.csv", index=False)
    pd.DataFrame(failures, columns=["analysis_set", "repeat", "reason"]).to_csv(
        ROOT / "results" / "fixed_taxon_panel_nested_budget_failures.csv", index=False
    )
    with open(ROOT / "results" / "fixed_taxon_panel_nested_budget_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "repeats": REPEATS,
            "test_taxon_fraction": TEST_TAXON_FRACTION,
            "budget_fractions_of_eligible_calibration_pairs": BUDGET_FRACTIONS,
            "test_taxa_and_test_pairs_fixed_across_budgets_within_repeat": True,
            "calibration_sets_nested_within_repeat": True,
            "calibration_pair_sampling_uses_outcomes": False,
            "partition_outcomes_used_only_to_require_both_test_classes": True,
            "model": "fixed L2 logistic regression, C=0.01, balanced class weights",
            "combined_features": FEATURE_SETS["combined"],
            "unsupervised_methods": list(METHODS),
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
