#!/usr/bin/env python3
"""Evaluate nested measurement budgets with five held-out pair folds."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from evaluate_direct_pair_features import FEATURE_SETS, load_analysis, score_model
from evaluate_within_system_calibration import ROOT, SEED


SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "friedman_microcosm_2017",
    "schafer_phyllosphere_2022",
]
BUDGET_FRACTIONS = [0.20, 0.40, 0.60, 0.80]
REFERENCE_METHODS = {
    "butyrate_assembly_2021": ("OneNet", "onenet_score"),
    "carlstrom_phyllosphere_2019": ("OneNet", "onenet_score"),
    "friedman_microcosm_2017": ("OneNet", "onenet_score"),
    "schafer_phyllosphere_2022": ("OneNet", "onenet_score"),
}


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    first = frame["taxon_1"].astype(str)
    second = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(first, second)
    frame["taxon_2"] = np.maximum(first, second)
    return frame


def load_system(system: str) -> pd.DataFrame:
    frame = load_analysis(system)
    onenet = canonicalize(pd.read_csv(
        ROOT / "analysis_data" / f"{system}_onenet_pair_scores.csv"
    ))
    frame = frame.merge(
        onenet[["dataset", "taxon_1", "taxon_2", "onenet_score"]],
        on=["dataset", "taxon_1", "taxon_2"],
        how="left",
        validate="many_to_one",
    )
    if frame["onenet_score"].isna().any():
        raise ValueError(f"{system}: missing OneNet scores after alignment")
    return frame


def run_system(system: str) -> list[dict]:
    frame = load_system(system)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    method_name, method_column = REFERENCE_METHODS[system]
    rows: list[dict] = []

    for fold, (pool_index, test_index) in enumerate(
        splitter.split(groups["pair_id"], groups["group_label"]), start=1
    ):
        pool_ids = groups.iloc[pool_index]["pair_id"].to_numpy()
        test_ids = set(groups.iloc[test_index]["pair_id"])
        test = frame[frame["pair_id"].isin(test_ids)].copy()
        y_test = test["interaction_label"].to_numpy()
        if np.unique(y_test).size < 2:
            raise RuntimeError(f"{system}, fold {fold}: test fold has one outcome class")

        reference_score = test[method_column].to_numpy()
        rows.append({
            "analysis_set": system,
            "fold": fold,
            "budget_fraction": 0.0,
            "budget_pairs": 0,
            "test_pairs": len(test_ids),
            "test_outcomes": len(test),
            "test_positives": int(y_test.sum()),
            "method": method_name,
            "auprc": average_precision_score(y_test, reference_score),
            "auroc": roc_auc_score(y_test, reference_score),
            "model_fitted": False,
        })

        rng = np.random.default_rng(SEED + 100_000 * SYSTEMS.index(system) + fold)
        ordered_pool = rng.permutation(pool_ids)
        total_pairs = len(groups)
        for fraction in BUDGET_FRACTIONS:
            requested = int(np.ceil(fraction * total_pairs))
            budget = min(requested, len(ordered_pool))
            training_ids = set(ordered_pool[:budget])
            training = frame[frame["pair_id"].isin(training_ids)].copy()
            _, score, fitted = score_model(training, test, FEATURE_SETS["combined"])
            rows.append({
                "analysis_set": system,
                "fold": fold,
                "budget_fraction": fraction,
                "budget_pairs": budget,
                "test_pairs": len(test_ids),
                "test_outcomes": len(test),
                "test_positives": int(y_test.sum()),
                "method": "Supervised combined",
                "auprc": average_precision_score(y_test, score),
                "auroc": roc_auc_score(y_test, score),
                "model_fitted": fitted,
            })
    return rows


def main() -> None:
    rows = []
    for system in SYSTEMS:
        rows.extend(run_system(system))
    results = pd.DataFrame(rows)
    summary = results.groupby(
        ["analysis_set", "method", "budget_fraction"], as_index=False
    ).agg(
        folds=("fold", "size"),
        minimum_budget_pairs=("budget_pairs", "min"),
        maximum_budget_pairs=("budget_pairs", "max"),
        mean_test_pairs=("test_pairs", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
        model_fit_rate=("model_fitted", "mean"),
    )
    results.to_csv(ROOT / "results" / "pair_five_fold_cv.csv", index=False)
    summary.to_csv(ROOT / "results" / "pair_five_fold_cv_summary.csv", index=False)
    with open(ROOT / "results" / "pair_five_fold_cv_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "folds": 5,
            "split_unit": "undirected pair identifier",
            "friedman_directions_kept_in_same_fold": True,
            "fold_assignment": "stratified by maximum pair outcome",
            "budget_fractions_of_all_tested_pair_identifiers": BUDGET_FRACTIONS,
            "training_subsets_nested_within_fold": True,
            "training_subset_order_uses_outcomes": False,
            "model": "fixed L2 logistic regression, C=0.01, balanced class weights",
            "combined_features": FEATURE_SETS["combined"],
            "reference_methods": {
                system: method for system, (method, _) in REFERENCE_METHODS.items()
            },
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
