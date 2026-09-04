#!/usr/bin/env python3
"""Compare taxon panel and distributed pair calibration inside one system."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

from evaluate_within_system_calibration import (
    FEATURES,
    ROOT,
    SEED,
    UNSUPERVISED,
    load_system,
    make_model,
)


SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
TAXON_FRACTIONS = [0.20, 0.40, 0.60]
REPEATS = 200


def pair_table(frame: pd.DataFrame) -> pd.DataFrame:
    """One row per pair for sampling; directional rows remain together."""
    return frame[["pair_id", "taxon_1", "taxon_2"]].drop_duplicates("pair_id")


def predict(calibration: pd.DataFrame, test: pd.DataFrame, seed: int):
    """Fit when two classes are observed; otherwise return the observed rate."""
    if calibration["interaction_label"].nunique() < 2:
        probability = float(calibration["interaction_label"].mean())
        return np.repeat(probability, len(test)), False
    model = make_model(0.01, True)
    model.fit(calibration[FEATURES], calibration["interaction_label"])
    predictions = model.predict_proba(test[FEATURES])[:, 1]
    return predictions, True


def evaluate_split(
    frame: pd.DataFrame,
    calibration_pairs: set[str],
    test_pairs: set[str],
    design: str,
    system: str,
    fraction: float,
    repeat: int,
    selected_taxa: int,
    seed: int,
) -> dict | None:
    calibration = frame[frame["pair_id"].isin(calibration_pairs)].copy()
    test = frame[frame["pair_id"].isin(test_pairs)].copy()
    if calibration.empty or test.empty or test["interaction_label"].nunique() < 2:
        return None
    scores, fitted = predict(calibration, test, seed)
    y = test["interaction_label"].to_numpy()
    row = {
        "analysis_set": system,
        "design": design,
        "taxon_fraction": fraction,
        "repeat": repeat,
        "selected_taxa": selected_taxa,
        "calibration_pairs": len(calibration_pairs),
        "calibration_outcomes": len(calibration),
        "calibration_positives": int(calibration["interaction_label"].sum()),
        "test_pairs": len(test_pairs),
        "test_outcomes": len(test),
        "test_positives": int(y.sum()),
        "test_prevalence": float(y.mean()),
        "model_fitted": fitted,
        "supervised_auprc": average_precision_score(y, scores),
        "mean_rank_auprc": average_precision_score(y, test["mean_score_rank"]),
    }
    row["supervised_minus_mean_rank"] = row["supervised_auprc"] - row["mean_rank_auprc"]
    for method, column in UNSUPERVISED.items():
        key = method.lower().replace(" ", "_").replace("-", "_")
        row[f"{key}_auprc"] = average_precision_score(y, test[column])
    return row


def run_system(system: str) -> tuple[list[dict], list[dict]]:
    frame = load_system(system)
    pairs = pair_table(frame)
    taxa = sorted(set(pairs["taxon_1"]) | set(pairs["taxon_2"]))
    rows, failures = [], []
    for fraction in TAXON_FRACTIONS:
        panel_size = max(4, int(round(len(taxa) * fraction)))
        panel_size = min(panel_size, len(taxa) - 3)
        for repeat in range(REPEATS):
            seed = SEED + 1_000_000 * SYSTEMS.index(system) + int(fraction * 10_000) + repeat
            rng = np.random.default_rng(seed)
            panel = set(rng.choice(taxa, size=panel_size, replace=False))
            panel_pairs = set(pairs.loc[
                pairs["taxon_1"].isin(panel) & pairs["taxon_2"].isin(panel), "pair_id"
            ])
            unseen_pairs = set(pairs.loc[
                ~pairs["taxon_1"].isin(panel) & ~pairs["taxon_2"].isin(panel), "pair_id"
            ])
            if not panel_pairs or not unseen_pairs:
                failures.append({"analysis_set": system, "design": "taxon_panel",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "empty calibration or test set"})
                continue
            panel_row = evaluate_split(
                frame, panel_pairs, unseen_pairs, "taxon_panel", system,
                fraction, repeat, panel_size, seed,
            )
            if panel_row is None:
                failures.append({"analysis_set": system, "design": "taxon_panel",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "test set lacks both classes"})
            else:
                rows.append(panel_row)

            budget = len(panel_pairs)
            if budget >= len(pairs):
                failures.append({"analysis_set": system, "design": "distributed_pairs",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "pair budget leaves no test set"})
                continue
            chosen = set(rng.choice(pairs["pair_id"].to_numpy(), size=budget, replace=False))
            remaining = set(pairs["pair_id"]) - chosen
            distributed_row = evaluate_split(
                frame, chosen, remaining, "distributed_pairs", system,
                fraction, repeat, panel_size, seed + 500_000,
            )
            if distributed_row is None:
                failures.append({"analysis_set": system, "design": "distributed_pairs",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "test set lacks both classes"})
            else:
                rows.append(distributed_row)
    return rows, failures


def summarize(results: pd.DataFrame) -> pd.DataFrame:
    summary = results.groupby(
        ["analysis_set", "design", "taxon_fraction"], as_index=False
    ).agg(
        valid_repeats=("repeat", "size"),
        mean_selected_taxa=("selected_taxa", "mean"),
        mean_calibration_pairs=("calibration_pairs", "mean"),
        calibration_positive_rate=("calibration_positives", lambda x: float((x > 0).mean())),
        model_fit_rate=("model_fitted", "mean"),
        mean_test_prevalence=("test_prevalence", "mean"),
        mean_rank_auprc=("mean_rank_auprc", "mean"),
        supervised_auprc=("supervised_auprc", "mean"),
        auprc_change=("supervised_minus_mean_rank", "mean"),
        supervised_win_rate=("supervised_minus_mean_rank", lambda x: float((x > 0).mean())),
    )
    return summary


def main() -> None:
    rows, failures = [], []
    for system in SYSTEMS:
        system_rows, system_failures = run_system(system)
        rows.extend(system_rows)
        failures.extend(system_failures)
    results = pd.DataFrame(rows)
    summary = summarize(results)
    results.to_csv(ROOT / "results" / "partial_system_calibration_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "partial_system_calibration_summary.csv", index=False)
    pd.DataFrame(failures).to_csv(
        ROOT / "results" / "partial_system_calibration_failures.csv", index=False
    )
    metadata = {
        "seed": SEED,
        "repeats": REPEATS,
        "taxon_fractions": TAXON_FRACTIONS,
        "features": FEATURES,
        "taxon_panel_training": "all scored tested pairs with both taxa in the sampled panel",
        "taxon_panel_test": "all scored tested pairs with both taxa outside the sampled panel",
        "distributed_training": "uniform random pairs across the scored tested universe",
        "distributed_test": "all remaining scored tested pairs",
        "label_blind_sampling": True,
        "distributed_pair_budget": "matched to the realized taxon panel pair count in each repeat",
        "one_class_training": "constant probability equal to the revealed positive fraction",
        "supervised_model": "L2 logistic regression with C=0.01 and balanced class weights, fixed before evaluation",
        "primary_baseline": "fixed mean percentile rank of PLNNetwork, Poisson GLMNet, and SparCC",
    }
    with open(ROOT / "results" / "partial_system_calibration_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    print(summary.to_string(index=False))
    print("\nFailures:", len(failures))


if __name__ == "__main__":
    main()
