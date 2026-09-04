#!/usr/bin/env python3
"""Add OneNet consensus scores to the fixed test benchmark."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from evaluate_fixed_test_nested_budgets import SYSTEMS
from evaluate_within_system_calibration import ROOT, canonicalize, load_system


def main() -> None:
    repeat_file = ROOT / "results" / "fixed_test_nested_budget_repeats.csv"
    repeats = pd.read_csv(repeat_file)
    template = repeats[repeats["method"] == "PLNNetwork"].copy()
    rows: list[dict] = []

    for system in SYSTEMS:
        frame = load_system(system)
        score_file = ROOT / "analysis_data" / f"{system}_onenet_pair_scores.csv"
        scores = canonicalize(pd.read_csv(score_file))
        merged = frame.merge(
            scores[["taxon_1", "taxon_2", "onenet_score"]],
            on=["taxon_1", "taxon_2"], how="left", validate="many_to_one",
        )
        if merged["onenet_score"].isna().any():
            missing = merged.loc[merged["onenet_score"].isna(), "pair_id"].unique()
            raise RuntimeError(f"{system}: {len(missing)} tested pairs lack OneNet scores")

        system_template = template[template["analysis_set"] == system]
        original = repeats[
            (repeats["analysis_set"] == system) & (repeats["method"] == "PLNNetwork")
        ]
        # Recover each fixed test set from the deterministic split used by the main analysis.
        from evaluate_fixed_test_nested_budgets import draw_fixed_test, TEST_FRACTION
        pair_ids = merged["pair_id"].drop_duplicates().to_numpy()
        test_size = max(2, int(np.ceil(TEST_FRACTION * len(pair_ids))))
        for repeat in sorted(system_template["repeat"].unique()):
            rng = np.random.default_rng(20260821 + 1_000_000 * SYSTEMS.index(system) + repeat)
            test_ids = draw_fixed_test(merged, pair_ids, test_size, rng)
            test = merged[merged["pair_id"].isin(test_ids)]
            y = test["interaction_label"].to_numpy()
            rows.append({
                "analysis_set": system,
                "repeat": repeat,
                "budget_fraction": 0.0,
                "budget_pairs": 0,
                "test_pairs": len(test_ids),
                "test_positives": int(y.sum()),
                "method": "OneNet (six estimators)",
                "auprc": average_precision_score(y, test["onenet_score"]),
                "auroc": roc_auc_score(y, test["onenet_score"]),
                "model_fitted": False,
            })

    updated = pd.concat([repeats[repeats["method"] != "OneNet"], pd.DataFrame(rows)], ignore_index=True)
    updated.to_csv(repeat_file, index=False)
    summary = updated.groupby(
        ["analysis_set", "method", "budget_fraction"], as_index=False
    ).agg(
        valid_repeats=("repeat", "size"), budget_pairs=("budget_pairs", "mean"),
        test_pairs=("test_pairs", "mean"), test_positives=("test_positives", "mean"),
        model_fit_rate=("model_fitted", "mean"), mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"), mean_auroc=("auroc", "mean"), sd_auroc=("auroc", "std"),
    )
    summary.to_csv(ROOT / "results" / "fixed_test_nested_budget_summary.csv", index=False)
    onenet = summary[summary["method"] == "OneNet (six estimators)"]
    onenet.to_csv(ROOT / "results" / "onenet_fixed_test_summary.csv", index=False)
    with open(ROOT / "results" / "onenet_fixed_test_metadata.json", "w") as handle:
        json.dump({
            "score": "mean edge selection frequency across six OneNet estimators",
            "mean_stability": 0.8,
            "resamples": 30,
            "estimators": ["PLNnetwork", "SpiecEasi", "gCoda", "EMtree", "Magma", "SPRING"],
            "omitted_estimator": "ZiLN; constant columns occurred in sparse resamples",
            "test_pairs": "same deterministic fixed test sets as the primary analysis",
            "experimental_labels_used_in_onenet": False,
        }, handle, indent=2)
    print(onenet.to_string(index=False))


if __name__ == "__main__":
    main()
