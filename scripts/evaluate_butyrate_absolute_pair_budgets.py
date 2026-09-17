#!/usr/bin/env python3
"""Repeated five-fold evaluation of absolute interaction-pair budgets in Butyrate."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from lens import ROOT, SEED, SYSTEMS, SYSTEM_SEED_INDEX, REFERENCE_METHODS, load_system, select_regularization, fit_score



SYSTEM = "butyrate_assembly_2021"
REPEATS = 20
FIXED_BUDGETS = [10, 20, 50]


def main() -> None:
    frame = load_system(SYSTEM)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    _, onenet_column = REFERENCE_METHODS[SYSTEM]
    system_index = SYSTEM_SEED_INDEX[SYSTEM]
    rows: list[dict] = []

    for repeat in range(1, REPEATS + 1):
        outer_seed = SEED + 10_000_000 * repeat
        splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=outer_seed)
        for fold, (pool_index, test_index) in enumerate(
            splitter.split(groups["pair_id"], groups["group_label"]), start=1
        ):
            pool_ids = groups.iloc[pool_index]["pair_id"].to_numpy()
            test_ids = set(groups.iloc[test_index]["pair_id"])
            test = frame[frame["pair_id"].isin(test_ids)].copy()
            y_test = test["interaction_label"].to_numpy(int)
            order_seed = outer_seed + 100_000 * system_index + fold
            rng = np.random.default_rng(order_seed)
            ordered_pool = rng.permutation(pool_ids)
            budgets = FIXED_BUDGETS + [len(ordered_pool)]

            for budget_index, budget in enumerate(budgets):
                training_ids = set(ordered_pool[:budget])
                training = frame[frame["pair_id"].isin(training_ids)].copy()
                tuning_seed = order_seed + 1_000 * budget_index
                candidate, _ = select_regularization(training, tuning_seed)
                score = fit_score(candidate, training, test, tuning_seed)
                onenet_score = test[onenet_column].to_numpy(float)
                supervised_auprc = average_precision_score(y_test, score)
                onenet_auprc = average_precision_score(y_test, onenet_score)
                rows.append({
                    "analysis_set": SYSTEM,
                    "repeat": repeat,
                    "fold": fold,
                    "budget_label": "maximum" if budget == len(ordered_pool) else str(budget),
                    "budget_pairs": budget,
                    "training_positives": int(training["interaction_label"].sum()),
                    "training_neutrals": int((training["interaction_label"] == 0).sum()),
                    "test_pairs": len(test_ids),
                    "test_positives": int(y_test.sum()),
                    "selected_setting": candidate.setting,
                    "supervised_auprc": supervised_auprc,
                    "onenet_auprc": onenet_auprc,
                    "auprc_difference": supervised_auprc - onenet_auprc,
                    "supervised_auroc": roc_auc_score(y_test, score),
                    "onenet_auroc": roc_auc_score(y_test, onenet_score),
                })

    results = pd.DataFrame(rows)
    repeat_means = results.groupby(
        ["repeat", "budget_label"], as_index=False
    ).agg(
        budget_pairs=("budget_pairs", "mean"),
        supervised_auprc=("supervised_auprc", "mean"),
        onenet_auprc=("onenet_auprc", "mean"),
        auprc_difference=("auprc_difference", "mean"),
    )
    summary = repeat_means.groupby("budget_label", as_index=False).agg(
        budget_pairs=("budget_pairs", "mean"),
        repeats=("repeat", "size"),
        mean_supervised_auprc=("supervised_auprc", "mean"),
        sd_supervised_auprc=("supervised_auprc", "std"),
        se_supervised_auprc=("supervised_auprc", lambda x: x.std() / np.sqrt(len(x))),
        mean_onenet_auprc=("onenet_auprc", "mean"),
        sd_onenet_auprc=("onenet_auprc", "std"),
        mean_auprc_difference=("auprc_difference", "mean"),
        sd_auprc_difference=("auprc_difference", "std"),
        se_auprc_difference=("auprc_difference", lambda x: x.std() / np.sqrt(len(x))),
        repeat_win_frequency=("auprc_difference", lambda x: np.mean(x > 0)),
    )
    order = {"10": 10, "20": 20, "50": 50, "maximum": 10_000}
    summary["order"] = summary["budget_label"].map(order)
    summary = summary.sort_values("order").drop(columns="order")

    results.to_csv(
        ROOT / "results" / "butyrate_absolute_pair_budgets_folds.csv", index=False
    )
    repeat_means.to_csv(
        ROOT / "results" / "butyrate_absolute_pair_budgets_repeats.csv", index=False
    )
    summary.to_csv(
        ROOT / "results" / "butyrate_absolute_pair_budgets_summary.csv", index=False
    )
    with open(
        ROOT / "results" / "butyrate_absolute_pair_budgets_metadata.json", "w"
    ) as handle:
        json.dump({
            "seed": SEED,
            "repeats": REPEATS,
            "folds_per_repeat": 5,
            "outer_split": "stratified by interaction outcome",
            "budget_pair_selection": "nested random order without outcomes",
            "fixed_budgets": FIXED_BUDGETS,
            "maximum_budget": "all pair identifiers outside each test fold",
            "model": "inner-tuned centered ridge logistic regression",
            "inner_validation_fallback": "mean percentile across PLNNetwork, Poisson GLMNet, and SparCC when either training outcome class has fewer than three pairs",
            "reference": "OneNet on the identical test pairs",
            "se_unit": "mean of five folds within each repeated partition",
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
