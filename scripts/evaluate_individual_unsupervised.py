#!/usr/bin/env python3
"""Evaluate each unsupervised score on the partial calibration test pairs."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from evaluate_within_system_calibration import ROOT, SEED, UNSUPERVISED, load_system


METHODS = {
    **UNSUPERVISED,
    "SPIEC-EASI": "spieceasi_score_percentile",
}


SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
TAXON_FRACTIONS = [0.20, 0.40, 0.60]
REPEATS = 200


def evaluate_test(frame: pd.DataFrame, test_pairs: set[str], system: str,
                  design: str, fraction: float, repeat: int,
                  calibration_pairs: int) -> list[dict]:
    test = frame[frame.pair_id.isin(test_pairs)]
    if test.empty or test.interaction_label.nunique() < 2:
        return []
    y = test.interaction_label.to_numpy()
    rows = []
    for method, column in METHODS.items():
        score = test[column].to_numpy()
        rows.append({
            "analysis_set": system,
            "design": design,
            "taxon_fraction": fraction,
            "repeat": repeat,
            "calibration_pairs": calibration_pairs,
            "test_pairs": test.pair_id.nunique(),
            "test_positives": int(y.sum()),
            "method": method,
            "auprc": average_precision_score(y, score),
            "auroc": roc_auc_score(y, score),
        })
    return rows


def run_system(system: str) -> tuple[list[dict], list[dict]]:
    frame = load_system(system)
    pairs = frame[["pair_id", "taxon_1", "taxon_2"]].drop_duplicates("pair_id")
    taxa = sorted(set(pairs.taxon_1) | set(pairs.taxon_2))
    rows, failures = [], []
    for fraction in TAXON_FRACTIONS:
        panel_size = min(max(4, int(round(len(taxa) * fraction))), len(taxa) - 3)
        for repeat in range(REPEATS):
            seed = SEED + 1_000_000 * SYSTEMS.index(system) + int(fraction * 10_000) + repeat
            rng = np.random.default_rng(seed)
            panel = set(rng.choice(taxa, size=panel_size, replace=False))
            panel_pairs = set(pairs.loc[
                pairs.taxon_1.isin(panel) & pairs.taxon_2.isin(panel), "pair_id"
            ])
            unseen_pairs = set(pairs.loc[
                ~pairs.taxon_1.isin(panel) & ~pairs.taxon_2.isin(panel), "pair_id"
            ])
            if not panel_pairs or not unseen_pairs:
                failures.append({"analysis_set": system, "design": "taxon_panel",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "empty calibration or test set"})
                continue
            panel_rows = evaluate_test(
                frame, unseen_pairs, system, "taxon_panel", fraction, repeat,
                len(panel_pairs),
            )
            if not panel_rows:
                failures.append({"analysis_set": system, "design": "taxon_panel",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "test set lacks both classes"})
            rows.extend(panel_rows)
            budget = len(panel_pairs)
            if budget >= len(pairs):
                continue
            chosen = set(rng.choice(pairs.pair_id.to_numpy(), size=budget, replace=False))
            distributed_rows = evaluate_test(
                frame, set(pairs.pair_id) - chosen, system, "distributed_pairs",
                fraction, repeat, budget,
            )
            if not distributed_rows:
                failures.append({"analysis_set": system, "design": "distributed_pairs",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "test set lacks both classes"})
            rows.extend(distributed_rows)
    return rows, failures


def main() -> None:
    rows, failures = [], []
    for system in SYSTEMS:
        system_rows, system_failures = run_system(system)
        rows.extend(system_rows)
        failures.extend(system_failures)
    results = pd.DataFrame(rows)
    summary = results.groupby(
        ["analysis_set", "design", "taxon_fraction", "method"], as_index=False
    ).agg(
        valid_repeats=("repeat", "size"),
        mean_calibration_pairs=("calibration_pairs", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
    )
    results.to_csv(ROOT / "results" / "individual_unsupervised_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "individual_unsupervised_summary.csv", index=False)
    pd.DataFrame(failures, columns=[
        "analysis_set", "design", "taxon_fraction", "repeat", "reason"
    ]).to_csv(ROOT / "results" / "individual_unsupervised_failures.csv", index=False)
    with open(ROOT / "results" / "individual_unsupervised_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "repeats": REPEATS,
            "taxon_fractions": TAXON_FRACTIONS,
            "methods": list(METHODS),
            "interaction_labels_used_for_method_fitting_or_selection": False,
            "test_pairs": "identical to the supervised partial calibration evaluation",
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
