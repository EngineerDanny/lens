#!/usr/bin/env python3
"""Export held-out predictions for the final 80% centred ridge analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

from evaluate_optimized_sparse_logistic import select_regularization
from evaluate_pair_five_fold_cv import REFERENCE_METHODS, SYSTEMS, load_system
from evaluate_within_system_calibration import ROOT, SEED
from optimize_supervised_model import fit_score


ANALYSIS_SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
]
BUDGET_FRACTION = 0.80


def run_system(system: str) -> list[dict]:
    frame = load_system(system)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    _, onenet_column = REFERENCE_METHODS[system]
    system_index = SYSTEMS.index(system)
    rows: list[dict] = []

    for fold, (pool_index, test_index) in enumerate(
        splitter.split(groups["pair_id"], groups["group_label"]), start=1
    ):
        pool_ids = groups.iloc[pool_index]["pair_id"].to_numpy()
        test_ids = set(groups.iloc[test_index]["pair_id"])
        test = frame[frame["pair_id"].isin(test_ids)].copy()
        rng = np.random.default_rng(SEED + 100_000 * system_index + fold)
        ordered_pool = rng.permutation(pool_ids)
        budget = min(int(np.ceil(BUDGET_FRACTION * len(groups))), len(ordered_pool))
        training_ids = set(ordered_pool[:budget])
        training = frame[frame["pair_id"].isin(training_ids)].copy()
        tuning_seed = (
            SEED + 1_000_000 * system_index + 10_000 * fold
            + int(100 * BUDGET_FRACTION)
        )
        candidate, _ = select_regularization(training, tuning_seed)
        supervised_score = fit_score(candidate, training, test, tuning_seed)

        for position, (_, row) in enumerate(test.iterrows()):
            rows.append({
                "analysis_set": system,
                "fold": fold,
                "pair_id": row["pair_id"],
                "taxon_1": row["taxon_1"],
                "taxon_2": row["taxon_2"],
                "interaction_label": int(row["interaction_label"]),
                "budget_fraction": BUDGET_FRACTION,
                "budget_pairs": budget,
                "selected_setting": candidate.setting,
                "supervised_score": float(supervised_score[position]),
                "onenet_score": float(row[onenet_column]),
            })
    return rows


def main() -> None:
    rows = [row for system in ANALYSIS_SYSTEMS for row in run_system(system)]
    output = pd.DataFrame(rows)
    output.to_csv(
        ROOT / "results" / "final_sparse_pr_predictions_80pct.csv", index=False
    )
    print(output.groupby("analysis_set").agg(
        tested_pairs=("pair_id", "nunique"),
        positives=("interaction_label", "sum"),
        folds=("fold", "nunique"),
    ).to_string())


if __name__ == "__main__":
    main()
