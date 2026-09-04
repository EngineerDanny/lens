#!/usr/bin/env python3
"""Nested leave-one-system-out search for supervised score calibration."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260821
RELEVANT = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
]
SCORE3 = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
SCORE5 = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "spieceasi_score_percentile",
    "spring_score_percentile",
    "sparcc_score_percentile",
]
OBS = [
    "prevalence_min",
    "prevalence_max",
    "joint_prevalence",
    "presence_jaccard",
    "overall_zero_frequency",
    "n_samples",
    "n_taxa",
]


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    a = frame["taxon_1"].astype(str)
    b = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(a, b)
    frame["taxon_2"] = np.maximum(a, b)
    return frame


def load_data() -> pd.DataFrame:
    frames = [canonicalize(pd.read_csv(ROOT / "training_data" / "all_labeled_pairs.csv"))]
    for dataset in ["carlstrom_phyllosphere_2019", "schafer_phyllosphere_2022"]:
        truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / f"{dataset}_tested_pairs.csv"))
        features = canonicalize(pd.read_csv(ROOT / "analysis_data" / f"{dataset}_pair_features.csv"))
        merged = truth.merge(
            features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one"
        )
        merged = merged[
            (merged["tested_status"].astype(str).str.lower() != "ambiguous")
            & merged["interaction_label"].notna()
        ]
        frames.append(merged)
    data = pd.concat(frames, ignore_index=True, sort=False)
    data["interaction_label"] = data["interaction_label"].astype(int)
    for score in SCORE5:
        for observable in ["prevalence_min", "joint_prevalence", "presence_jaccard"]:
            data[f"{score}_x_{observable}"] = data[score] * data[observable]
    return data


def study_weights(data: pd.DataFrame, balance_classes: bool) -> np.ndarray:
    weights = np.zeros(len(data), dtype=float)
    for _, indices in data.groupby("dataset").groups.items():
        indices = np.asarray(list(indices), dtype=int)
        y = data.loc[indices, "interaction_label"].to_numpy()
        if balance_classes and len(np.unique(y)) == 2:
            for value in [0, 1]:
                mask = y == value
                weights[indices[mask]] = 0.5 / mask.sum()
        else:
            weights[indices] = 1 / len(indices)
    return weights / weights.mean()


def candidate_specs() -> list[dict]:
    interactions = [c for c in load_data().columns if "_x_" in c]
    feature_sets = {
        "score3": SCORE3,
        "score5": SCORE5,
        "score3_obs": SCORE3 + OBS,
        "score5_obs": SCORE5 + OBS,
        "score5_interactions": SCORE5 + OBS + interactions,
    }
    specs = []
    for feature_name, features in feature_sets.items():
        for c_value in [0.001, 0.01, 0.1, 1.0, 10.0]:
            for balanced in [False, True]:
                specs.append({
                    "name": f"logistic_{feature_name}_C{c_value}_balanced{balanced}",
                    "kind": "logistic", "features": features, "C": c_value,
                    "balanced": balanced,
                })
    for depth in [2, 4, 6]:
        for leaf in [10, 30]:
            specs.append({
                "name": f"extra_score5_obs_depth{depth}_leaf{leaf}",
                "kind": "extra", "features": SCORE5 + OBS, "depth": depth,
                "leaf": leaf, "balanced": False,
            })
    for leaves in [3, 7, 15]:
        for l2 in [0.1, 1.0, 10.0]:
            specs.append({
                "name": f"hist_score5_obs_leaves{leaves}_l2{l2}",
                "kind": "hist", "features": SCORE5 + OBS, "leaves": leaves,
                "l2": l2, "balanced": False,
            })
    return specs


def make_model(spec: dict):
    if spec["kind"] == "logistic":
        return Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("model", LogisticRegression(
                penalty="l2", solver="liblinear", C=spec["C"], max_iter=5000,
                random_state=SEED,
            )),
        ])
    if spec["kind"] == "extra":
        return Pipeline([
            ("impute", SimpleImputer(strategy="median")),
            ("model", ExtraTreesClassifier(
                n_estimators=300, max_depth=spec["depth"],
                min_samples_leaf=spec["leaf"], max_features="sqrt",
                random_state=SEED, n_jobs=-1,
            )),
        ])
    return Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("model", HistGradientBoostingClassifier(
            max_iter=150, max_leaf_nodes=spec["leaves"],
            l2_regularization=spec["l2"], learning_rate=0.05,
            random_state=SEED,
        )),
    ])


def fit_predict(train: pd.DataFrame, test: pd.DataFrame, spec: dict) -> np.ndarray:
    model = make_model(spec)
    weights = study_weights(train.reset_index(drop=True), spec["balanced"])
    model.fit(train[spec["features"]], train["interaction_label"], model__sample_weight=weights)
    return model.predict_proba(test[spec["features"]])[:, 1]


def inner_score(train: pd.DataFrame, spec: dict) -> tuple[float, float, int]:
    lifts = []
    for held_out in sorted(train["dataset"].unique()):
        validation = train[train["dataset"] == held_out]
        fitting = train[train["dataset"] != held_out]
        if validation["interaction_label"].nunique() < 2 or fitting["interaction_label"].nunique() < 2:
            continue
        score = fit_predict(fitting, validation, spec)
        y = validation["interaction_label"].to_numpy()
        lifts.append(average_precision_score(y, score) / y.mean())
    return float(np.mean(lifts)), float(np.median(lifts)), len(lifts)


def main() -> None:
    data = load_data()
    specs = candidate_specs()
    tuning_rows = []
    result_rows = []
    prediction_frames = []

    for held_out in RELEVANT:
        train = data[data["dataset"] != held_out].reset_index(drop=True)
        test = data[data["dataset"] == held_out].reset_index(drop=True)
        scored_specs = []
        for spec in specs:
            mean_lift, median_lift, folds = inner_score(train, spec)
            tuning_rows.append({
                "held_out": held_out, "candidate": spec["name"],
                "inner_mean_auprc_lift": mean_lift,
                "inner_median_auprc_lift": median_lift, "inner_folds": folds,
            })
            scored_specs.append((mean_lift, median_lift, spec["name"], spec))
        selected = max(scored_specs, key=lambda item: (item[0], item[1], item[2]))[3]
        score = fit_predict(train, test, selected)
        y = test["interaction_label"].to_numpy()
        result_rows.append({
            "analysis_set": held_out, "scope": "full", "selected_model": selected["name"],
            "tested_pairs": len(test), "positive_pairs": int(y.sum()),
            "random_auprc": float(y.mean()),
            "supervised_auprc": average_precision_score(y, score),
        })
        predictions = test[["dataset", "taxon_1", "taxon_2", "interaction_label"]].copy()
        predictions["supervised_score"] = score
        predictions["selected_model"] = selected["name"]
        prediction_frames.append(predictions)

        if held_out == "schafer_phyllosphere_2022":
            eligible = canonicalize(pd.read_csv(
                ROOT / "analysis_data" / f"{held_out}_conet_tested_pair_scores.csv"
            ))
            eligible_pairs = set(zip(eligible["taxon_1"], eligible["taxon_2"]))
            mask = np.array([
                (a, b) in eligible_pairs for a, b in zip(test["taxon_1"], test["taxon_2"])
            ])
            restricted_y = y[mask]
            result_rows.append({
                "analysis_set": held_out, "scope": "conet_universe",
                "selected_model": selected["name"], "tested_pairs": int(mask.sum()),
                "positive_pairs": int(restricted_y.sum()),
                "random_auprc": float(restricted_y.mean()),
                "supervised_auprc": average_precision_score(restricted_y, score[mask]),
            })

    results = pd.DataFrame(result_rows)
    tuning = pd.DataFrame(tuning_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    results.to_csv(ROOT / "results" / "nested_supervised_search_performance.csv", index=False)
    tuning.to_csv(ROOT / "results" / "nested_supervised_search_tuning.csv", index=False)
    predictions.to_csv(ROOT / "analysis_data" / "nested_supervised_search_predictions.csv", index=False)
    metadata = {
        "random_seed": SEED,
        "outer_validation": "leave one complete experimental system out",
        "inner_selection": "highest mean AUPRC lift across held-out training systems",
        "held_out_labels_used_for_selection": False,
        "candidate_count": len(specs),
    }
    with open(ROOT / "results" / "nested_supervised_search_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
