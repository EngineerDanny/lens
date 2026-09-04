#!/usr/bin/env python3
"""Minimal comparison of random and stratified taxon-covering budget sampling."""

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
    "schafer_phyllosphere_2022",
]
ORIGINAL_SYSTEM_INDEX = {
    "butyrate_assembly_2021": 0,
    "carlstrom_phyllosphere_2019": 1,
    "schafer_phyllosphere_2022": 3,
}
BUDGET_FRACTIONS = [0.20, 0.40, 0.60, 0.80]


def taxon_covering_order(groups: pd.DataFrame, seed: int) -> list[str]:
    """Return a nested, label-blind order favoring broad early taxon coverage."""
    rng = np.random.default_rng(seed)
    work = groups.copy()
    work["tie"] = rng.random(len(work))
    remaining = set(work.index)
    selected: list[str] = []
    covered: set[str] = set()
    while remaining:
        candidates = work.loc[list(remaining)]
        candidates = candidates.assign(
            new_taxa=candidates.apply(
                lambda row: int(row["taxon_1"] not in covered)
                + int(row["taxon_2"] not in covered), axis=1
            )
        ).sort_values(["new_taxa", "tie"], ascending=[False, True])
        chosen_index = int(candidates.index[0])
        chosen = work.loc[chosen_index]
        selected.append(str(chosen["pair_id"]))
        covered.update([str(chosen["taxon_1"]), str(chosen["taxon_2"])])
        remaining.remove(chosen_index)
    return selected


def run_system(system: str) -> list[dict]:
    frame = load_analysis(system)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max"),
        taxon_1=("taxon_1", "first"),
        taxon_2=("taxon_2", "first"),
    )
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    rows: list[dict] = []

    for fold, (pool_index, test_index) in enumerate(
        splitter.split(groups["pair_id"], groups["group_label"]), start=1
    ):
        pool = groups.iloc[pool_index].reset_index(drop=True)
        test_ids = set(groups.iloc[test_index]["pair_id"])
        test = frame[frame["pair_id"].isin(test_ids)].copy()
        y_test = test["interaction_label"].to_numpy()
        # Preserve the seeds used by the original four-system experiment.
        seed = SEED + 100_000 * ORIGINAL_SYSTEM_INDEX[system] + fold
        rng = np.random.default_rng(seed)
        orders = {
            "Old random": list(rng.permutation(pool["pair_id"].to_numpy())),
            "Taxon coverage": taxon_covering_order(pool, seed),
        }
        total_pairs = len(groups)
        for design, order in orders.items():
            for fraction in BUDGET_FRACTIONS:
                budget = min(int(np.ceil(fraction * total_pairs)), len(order))
                training_ids = set(order[:budget])
                training = frame[frame["pair_id"].isin(training_ids)].copy()
                _, score, fitted = score_model(training, test, FEATURE_SETS["combined"])
                rows.append({
                    "analysis_set": system,
                    "fold": fold,
                    "design": design,
                    "budget_fraction": fraction,
                    "budget_pairs": budget,
                    "training_positives": int(training["interaction_label"].sum()),
                    "training_positive_fraction": float(training["interaction_label"].mean()),
                    "training_taxa": len(set(training["taxon_1"]) | set(training["taxon_2"])),
                    "test_pairs": len(test_ids),
                    "test_positives": int(y_test.sum()),
                    "auprc": average_precision_score(y_test, score),
                    "auroc": roc_auc_score(y_test, score),
                    "model_fitted": fitted,
                })
    return rows


def main() -> None:
    rows = [row for system in SYSTEMS for row in run_system(system)]
    results = pd.DataFrame(rows)
    summary = results.groupby(
        ["analysis_set", "design", "budget_fraction"], as_index=False
    ).agg(
        folds=("fold", "size"),
        budget_pairs=("budget_pairs", "median"),
        mean_training_positives=("training_positives", "mean"),
        mean_training_taxa=("training_taxa", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
    )
    comparison = summary.pivot(
        index=["analysis_set", "budget_fraction"], columns="design", values="mean_auprc"
    ).reset_index()
    comparison["auprc_change"] = (
        comparison["Taxon coverage"] - comparison["Old random"]
    )
    results.to_csv(ROOT / "results" / "budget_sampling_design_comparison_folds.csv", index=False)
    summary.to_csv(ROOT / "results" / "budget_sampling_design_comparison_summary.csv", index=False)
    comparison.to_csv(ROOT / "results" / "budget_sampling_design_comparison.csv", index=False)
    with open(ROOT / "results" / "budget_sampling_design_comparison_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "outer_folds": 5,
            "outer_fold_assignment": "stratified by pair outcome",
            "test_folds_identical_between_designs": True,
            "model_changed": False,
            "new_budget_sampling": "nested label-blind order with early taxon coverage",
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
