#!/usr/bin/env python3
"""Apply the frozen three-score calibration to the Friedman benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
DATASET = "friedman_microcosm_2017"
SEED = 20260821
PREDICTORS = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
METHODS = ["pln", "glmnet", "spieceasi", "spring", "sparcc"]


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    ordered = np.sort(frame[["taxon_1", "taxon_2"]].astype(str).to_numpy(), axis=1)
    frame[["taxon_1", "taxon_2"]] = ordered
    return frame


def development_data() -> pd.DataFrame:
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
    return pd.concat(frames, ignore_index=True, sort=False)


def study_weights(data: pd.DataFrame) -> np.ndarray:
    counts = data["dataset"].value_counts()
    weights = data["dataset"].map(1 / counts).to_numpy(float)
    return weights / weights.mean()


def clustered_ap_interval(data: pd.DataFrame, score_column: str) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    groups = [group for _, group in data.groupby(["taxon_1", "taxon_2"], sort=True)]
    values = []
    for _ in range(2000):
        sampled = [groups[index] for index in rng.integers(0, len(groups), len(groups))]
        bootstrap = pd.concat(sampled, ignore_index=True)
        y = bootstrap["interaction_label"].astype(int)
        if y.nunique() == 2:
            values.append(average_precision_score(y, bootstrap[score_column].astype(float)))
    return float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def main() -> None:
    train = development_data()
    pair_features = canonicalize(
        pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_pair_features.csv")
    )

    model = Pipeline([
        ("impute", SimpleImputer(strategy="median")),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2", solver="liblinear", C=1.0, max_iter=5000,
            random_state=SEED,
        )),
    ])
    model.fit(
        train[PREDICTORS], train["interaction_label"].astype(int),
        model__sample_weight=study_weights(train),
    )
    pair_predictions = pair_features[["dataset", "taxon_1", "taxon_2"]].copy()
    pair_predictions["supervised_score"] = model.predict_proba(pair_features[PREDICTORS])[:, 1]
    for method in METHODS:
        pair_predictions[f"{method}_score"] = pair_features[f"{method}_score_percentile"]

    conet_path = ROOT / "analysis_data" / f"{DATASET}_conet_tested_pair_scores.csv"
    if conet_path.exists():
        conet = canonicalize(pd.read_csv(conet_path))
        pair_predictions = pair_predictions.merge(
            conet[["dataset", "taxon_1", "taxon_2", "conet_score"]],
            on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one",
        )

    # Persist predictions before opening the experimental outcome table.
    prediction_path = ROOT / "analysis_data" / f"{DATASET}_frozen_predictions.csv"
    pair_predictions.to_csv(prediction_path, index=False)

    truth = pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_tested_pairs_directional.csv")
    truth["focal_taxon"] = truth["taxon_1"]
    truth["partner_taxon"] = truth["taxon_2"]
    truth = canonicalize(truth)
    evaluated = truth.merge(
        pair_predictions, on=["dataset", "taxon_1", "taxon_2"], validate="many_to_one"
    )
    y = evaluated["interaction_label"].astype(int)
    rows = []
    evaluation_methods = METHODS + (["conet"] if "conet_score" in evaluated else []) + ["supervised"]
    for method in evaluation_methods:
        score_column = "supervised_score" if method == "supervised" else f"{method}_score"
        score = evaluated[score_column].astype(float)
        ci_low, ci_high = clustered_ap_interval(evaluated, score_column)
        rows.append({
            "method": method,
            "tested_directions": len(evaluated),
            "distinct_tested_pairs": evaluated[["taxon_1", "taxon_2"]].drop_duplicates().shape[0],
            "supported_effects": int(y.sum()),
            "tested_neutrals": int((y == 0).sum()),
            "random_auprc": float(y.mean()),
            "average_precision": average_precision_score(y, score),
            "cluster_bootstrap_ap_ci_low": ci_low,
            "cluster_bootstrap_ap_ci_high": ci_high,
            "auroc": roc_auc_score(y, score),
        })
    performance = pd.DataFrame(rows)
    performance.to_csv(ROOT / "results" / f"{DATASET}_confirmatory_performance.csv", index=False)
    evaluated.to_csv(ROOT / "analysis_data" / f"{DATASET}_evaluated_directional_predictions.csv", index=False)

    coefficients = pd.DataFrame({
        "predictor": PREDICTORS,
        "standardized_coefficient": model.named_steps["model"].coef_[0],
    })
    coefficients.to_csv(ROOT / "results" / f"{DATASET}_frozen_coefficients.csv", index=False)
    metadata = {
        "model": "L2 logistic regression",
        "C": 1.0,
        "class_weight": None,
        "training_system_weighting": "equal total weight per experimental system",
        "predictors": PREDICTORS,
        "development_systems": sorted(train["dataset"].unique().tolist()),
        "development_rows": len(train),
        "friedman_labels_used_for_fitting_or_selection": False,
        "prediction_file_written_before_truth_loaded": True,
        "random_seed": SEED,
    }
    with open(ROOT / "results" / f"{DATASET}_confirmatory_metadata.json", "w") as handle:
        json.dump(metadata, handle, indent=2)
    print(performance.to_string(index=False))
    print("\nFrozen coefficients")
    print(coefficients.to_string(index=False))


if __name__ == "__main__":
    main()
