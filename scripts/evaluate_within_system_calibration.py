#!/usr/bin/env python3
"""Evaluate calibration from a subset of labels in the target microbial system."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedGroupKFold, StratifiedShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260821
REPEATS = 200
BUDGETS = [0.20, 0.40, 0.60]
SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
FEATURES = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
UNSUPERVISED = {
    "PLNNetwork": "pln_score_percentile",
    "Poisson GLMNet": "glmnet_score_percentile",
    "SparCC": "sparcc_score_percentile",
    "Mean rank": "mean_score_rank",
}
C_VALUES = [0.001, 0.01, 0.1, 1.0, 10.0]


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    first = frame["taxon_1"].astype(str)
    second = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(first, second)
    frame["taxon_2"] = np.maximum(first, second)
    return frame


def load_system(system: str) -> pd.DataFrame:
    if system == "butyrate_assembly_2021":
        frame = pd.read_csv(ROOT / "training_data" / "all_labeled_pairs.csv")
        frame = frame[frame["dataset"] == system].copy()
    else:
        truth_name = (
            f"{system}_tested_pairs_directional.csv"
            if system == "friedman_microcosm_2017"
            else f"{system}_tested_pairs.csv"
        )
        truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / truth_name))
        features = canonicalize(pd.read_csv(
            ROOT / "analysis_data" / f"{system}_pair_features.csv"
        ))
        frame = truth.merge(
            features, on=["dataset", "taxon_1", "taxon_2"],
            how="inner", validate="many_to_one",
        )
    frame = canonicalize(frame)
    frame = frame[frame["interaction_label"].notna()].copy().reset_index(drop=True)
    frame["interaction_label"] = frame["interaction_label"].astype(int)
    frame["pair_id"] = frame["taxon_1"] + "||" + frame["taxon_2"]
    frame["mean_score_rank"] = frame[FEATURES].mean(axis=1)
    return frame


def make_model(c_value: float, balanced: bool) -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2", solver="liblinear", C=c_value,
            class_weight="balanced" if balanced else None,
            max_iter=5000, random_state=SEED,
        )),
    ])


def tune_model(calibration: pd.DataFrame, seed: int) -> tuple[float, bool]:
    y = calibration["interaction_label"].to_numpy()
    group_labels = calibration.groupby("pair_id")["interaction_label"].max()
    group_minority = min(int(group_labels.sum()), int((1 - group_labels).sum()))
    if group_minority < 2:
        return 0.01, True
    folds = min(5, group_minority)
    cv = StratifiedGroupKFold(n_splits=folds, shuffle=True, random_state=seed)
    candidates = []
    for c_value in C_VALUES:
        for balanced in [False, True]:
            score = np.full(len(calibration), np.nan)
            for fitting_index, validation_index in cv.split(
                    calibration[FEATURES], y, calibration["pair_id"]):
                model = make_model(c_value, balanced)
                model.fit(calibration.iloc[fitting_index][FEATURES], y[fitting_index])
                score[validation_index] = model.predict_proba(
                    calibration.iloc[validation_index][FEATURES]
                )[:, 1]
            if np.isnan(score).any():
                continue
            candidates.append((average_precision_score(y, score), -c_value,
                               not balanced, c_value, balanced))
    if not candidates:
        return 0.01, True
    chosen = max(candidates)
    return chosen[3], chosen[4]


def fit_supervised(calibration: pd.DataFrame, test: pd.DataFrame,
                   seed: int) -> tuple[np.ndarray, float, bool]:
    c_value, balanced = tune_model(calibration, seed)
    model = make_model(c_value, balanced)
    model.fit(calibration[FEATURES], calibration["interaction_label"])
    return model.predict_proba(test[FEATURES])[:, 1], c_value, balanced


def select_unsupervised(calibration: pd.DataFrame) -> str:
    y = calibration["interaction_label"]
    prevalence = y.mean()
    candidates = []
    for method, column in UNSUPERVISED.items():
        score = average_precision_score(y, calibration[column]) / prevalence
        candidates.append((score, method))
    return max(candidates)[1]


def pair_groups(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )


def valid_split(frame: pd.DataFrame, calibration_pairs: set[str]) -> bool:
    calibration = frame[frame["pair_id"].isin(calibration_pairs)]
    test = frame[~frame["pair_id"].isin(calibration_pairs)]
    return (calibration["interaction_label"].nunique() == 2
            and test["interaction_label"].nunique() == 2)


def draw_split(frame: pd.DataFrame, fraction: float, seed: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    groups = pair_groups(frame)
    for attempt in range(100):
        splitter = StratifiedShuffleSplit(
            n_splits=1, train_size=fraction, random_state=seed + attempt
        )
        calibration_index, _ = next(splitter.split(groups, groups["group_label"]))
        selected = set(groups.iloc[calibration_index]["pair_id"])
        if valid_split(frame, selected):
            return (frame[frame["pair_id"].isin(selected)].copy(),
                    frame[~frame["pair_id"].isin(selected)].copy())
    raise RuntimeError("Could not obtain a split containing both classes")


def evaluate_pair_holdout(frame: pd.DataFrame, system: str) -> list[dict]:
    rows = []
    for budget in BUDGETS:
        for repeat in range(REPEATS):
            split_seed = SEED + 100000 * SYSTEMS.index(system) + 1000 * int(budget * 100) + repeat
            calibration, test = draw_split(frame, budget, split_seed)
            supervised, c_value, balanced = fit_supervised(calibration, test, split_seed)
            selected_method = select_unsupervised(calibration)
            selected_column = UNSUPERVISED[selected_method]
            y = test["interaction_label"].to_numpy()
            prevalence = y.mean()
            supervised_ap = average_precision_score(y, supervised)
            selected_ap = average_precision_score(y, test[selected_column])
            mean_rank_ap = average_precision_score(y, test["mean_score_rank"])
            row = {
                "analysis_set": system, "label_budget": budget, "repeat": repeat,
                "calibration_pairs": calibration["pair_id"].nunique(),
                "test_pairs": test["pair_id"].nunique(),
                "calibration_outcomes": len(calibration), "test_outcomes": len(test),
                "test_positives": int(y.sum()), "test_prevalence": prevalence,
                "selected_unsupervised": selected_method,
                "selected_unsupervised_auprc": selected_ap,
                "mean_rank_auprc": mean_rank_ap,
                "supervised_auprc": supervised_ap,
                "supervised_minus_selected": supervised_ap - selected_ap,
                "supervised_minus_mean_rank": supervised_ap - mean_rank_ap,
                "selected_C": c_value, "balanced_classes": balanced,
            }
            for method, column in UNSUPERVISED.items():
                key = method.lower().replace(" ", "_").replace("-", "_")
                row[f"{key}_auprc"] = average_precision_score(y, test[column])
            rows.append(row)
    return rows


def evaluate_taxon_holdout(frame: pd.DataFrame, system: str) -> list[dict]:
    rows = []
    taxa = sorted(set(frame["taxon_1"]) | set(frame["taxon_2"]))
    for taxon_index, taxon in enumerate(taxa):
        test_mask = (frame["taxon_1"] == taxon) | (frame["taxon_2"] == taxon)
        calibration, test = frame[~test_mask].copy(), frame[test_mask].copy()
        if (calibration["interaction_label"].nunique() < 2
                or test["interaction_label"].nunique() < 2
                or min(test["interaction_label"].value_counts()) < 2):
            continue
        seed = SEED + 10000 * SYSTEMS.index(system) + taxon_index
        supervised, c_value, balanced = fit_supervised(calibration, test, seed)
        selected_method = select_unsupervised(calibration)
        selected_column = UNSUPERVISED[selected_method]
        y = test["interaction_label"]
        supervised_ap = average_precision_score(y, supervised)
        selected_ap = average_precision_score(y, test[selected_column])
        rows.append({
            "analysis_set": system, "held_out_taxon": taxon,
            "test_outcomes": len(test), "test_positives": int(y.sum()),
            "test_prevalence": y.mean(), "selected_unsupervised": selected_method,
            "selected_unsupervised_auprc": selected_ap,
            "supervised_auprc": supervised_ap,
            "supervised_minus_selected": supervised_ap - selected_ap,
            "selected_C": c_value, "balanced_classes": balanced,
        })
    return rows


def summarize_pair_results(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (system, budget), group in results.groupby(["analysis_set", "label_budget"], sort=False):
        difference = group["supervised_minus_selected"]
        rows.append({
            "analysis_set": system, "label_budget": budget,
            "repeats": len(group),
            "median_calibration_pairs": int(group["calibration_pairs"].median()),
            "mean_selected_unsupervised_auprc": group["selected_unsupervised_auprc"].mean(),
            "mean_supervised_auprc": group["supervised_auprc"].mean(),
            "mean_auprc_difference": difference.mean(),
            "difference_ci_low": difference.quantile(0.025),
            "difference_ci_high": difference.quantile(0.975),
            "supervised_win_rate": (difference > 0).mean(),
        })
    return pd.DataFrame(rows)


def main() -> None:
    pair_rows, taxon_rows = [], []
    for system in SYSTEMS:
        frame = load_system(system)
        pair_rows.extend(evaluate_pair_holdout(frame, system))
        taxon_rows.extend(evaluate_taxon_holdout(frame, system))

    pair_results = pd.DataFrame(pair_rows)
    summary = summarize_pair_results(pair_results)
    taxon_results = pd.DataFrame(taxon_rows)
    pair_results.to_csv(ROOT / "results" / "within_system_calibration_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "within_system_calibration_summary.csv", index=False)
    taxon_results.to_csv(ROOT / "results" / "within_system_taxon_holdout.csv", index=False)
    with open(ROOT / "results" / "within_system_calibration_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED, "repeats": REPEATS, "label_budgets": BUDGETS,
            "features": FEATURES,
            "pair_split_unit": "undirected taxon pair; directional outcomes remain together",
            "model_selection": "C and class weighting selected by stratified cross validation using revealed labels only",
            "test_labels_used_for_model_selection": False,
        }, handle, indent=2)
    print(summary.to_string(index=False))
    print("\nTaxon holdout\n", taxon_results.groupby("analysis_set").agg(
        folds=("held_out_taxon", "size"),
        selected_auprc=("selected_unsupervised_auprc", "mean"),
        supervised_auprc=("supervised_auprc", "mean"),
        mean_difference=("supervised_minus_selected", "mean"),
    ).to_string())


if __name__ == "__main__":
    main()
