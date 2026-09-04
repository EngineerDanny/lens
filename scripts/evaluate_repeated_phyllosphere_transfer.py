#!/usr/bin/env python3
"""Transfer a score calibrator between two Arabidopsis phyllosphere studies."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260821
STUDIES = ["carlstrom_phyllosphere_2019", "schafer_phyllosphere_2022"]
FEATURES = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
METHODS = {
    "PLNNetwork": "pln_score_percentile",
    "Poisson GLMNet": "glmnet_score_percentile",
    "SparCC": "sparcc_score_percentile",
    "Mean rank": "mean_rank",
}
C_VALUES = [0.001, 0.01, 0.1, 1.0, 10.0]


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    first = frame["taxon_1"].astype(str)
    second = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(first, second)
    frame["taxon_2"] = np.maximum(first, second)
    frame["pair_id"] = frame["taxon_1"] + "||" + frame["taxon_2"]
    return frame


def load_study(study: str) -> pd.DataFrame:
    truth = canonicalize(pd.read_csv(
        ROOT / "cleaned_data" / f"{study}_tested_pairs.csv"
    ))
    scores = canonicalize(pd.read_csv(
        ROOT / "analysis_data" / f"{study}_pair_features.csv"
    ))
    frame = truth.merge(
        scores.drop(columns=["dataset"]),
        on=["taxon_1", "taxon_2", "pair_id"],
        how="inner",
        validate="one_to_one",
    )
    frame = frame[frame["interaction_label"].notna()].copy()
    frame["interaction_label"] = frame["interaction_label"].astype(int)
    frame["mean_rank"] = frame[FEATURES].mean(axis=1)
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


def tune(train: pd.DataFrame) -> tuple[float, bool, float]:
    y = train["interaction_label"].to_numpy()
    folds = min(5, int(np.bincount(y).min()))
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=SEED)
    candidates = []
    for c_value in C_VALUES:
        for balanced in (False, True):
            predictions = np.full(len(train), np.nan)
            for fit_index, validation_index in cv.split(train[FEATURES], y):
                model = make_model(c_value, balanced)
                model.fit(train.iloc[fit_index][FEATURES], y[fit_index])
                predictions[validation_index] = model.predict_proba(
                    train.iloc[validation_index][FEATURES]
                )[:, 1]
            ap = average_precision_score(y, predictions)
            candidates.append((ap, -c_value, not balanced, c_value, balanced))
    chosen = max(candidates)
    return chosen[3], chosen[4], chosen[0]


def main() -> None:
    frames = {study: load_study(study) for study in STUDIES}
    raw_truth = {
        study: canonicalize(pd.read_csv(
            ROOT / "cleaned_data" / f"{study}_tested_pairs.csv"
        )) for study in STUDIES
    }
    experimentally_common = set(raw_truth[STUDIES[0]]["pair_id"]) & set(
        raw_truth[STUDIES[1]]["pair_id"]
    )
    common_pairs = set(frames[STUDIES[0]]["pair_id"]) & set(frames[STUDIES[1]]["pair_id"])
    common = {
        study: frame[frame["pair_id"].isin(common_pairs)].copy().sort_values("pair_id")
        for study, frame in frames.items()
    }

    labels = common[STUDIES[0]][["pair_id", "interaction_label"]].merge(
        common[STUDIES[1]][["pair_id", "interaction_label"]],
        on="pair_id", suffixes=("_carlstrom", "_schafer"), validate="one_to_one",
    )
    agreement = pd.crosstab(
        labels["interaction_label_carlstrom"],
        labels["interaction_label_schafer"],
    ).rename_axis("carlstrom_label").rename_axis("schafer_label", axis=1)

    rows = []
    predictions = []
    for source, target in ((STUDIES[0], STUDIES[1]), (STUDIES[1], STUDIES[0])):
        train, test = common[source], common[target]
        c_value, balanced, training_cv_ap = tune(train)
        model = make_model(c_value, balanced)
        model.fit(train[FEATURES], train["interaction_label"])
        supervised = model.predict_proba(test[FEATURES])[:, 1]
        target_y = test["interaction_label"].to_numpy()
        source_y = train["interaction_label"].to_numpy()
        method_scores = {
            method: average_precision_score(target_y, test[column])
            for method, column in METHODS.items()
        }
        best_method = max(method_scores, key=method_scores.get)
        rows.append({
            "training_study": source,
            "test_study": target,
            "common_pairs": len(common_pairs),
            "training_positives": int(source_y.sum()),
            "test_positives": int(target_y.sum()),
            "test_prevalence": float(target_y.mean()),
            "supervised_calibration_auprc": average_precision_score(target_y, supervised),
            "best_unsupervised_method": best_method,
            "best_unsupervised_auprc": method_scores[best_method],
            "pln_auprc": method_scores["PLNNetwork"],
            "glmnet_auprc": method_scores["Poisson GLMNet"],
            "sparcc_auprc": method_scores["SparCC"],
            "mean_rank_auprc": method_scores["Mean rank"],
            "selected_C": c_value,
            "balanced_classes": balanced,
            "training_cv_auprc": training_cv_ap,
        })
        output = test[["pair_id", "taxon_1", "taxon_2", "interaction_label"]].copy()
        output.insert(0, "training_study", source)
        output.insert(1, "test_study", target)
        output["supervised_score"] = supervised
        for method, column in METHODS.items():
            output[method.lower().replace(" ", "_")] = test[column].to_numpy()
        predictions.append(output)

    result = pd.DataFrame(rows)
    result["supervised_minus_best_unsupervised"] = (
        result["supervised_calibration_auprc"] - result["best_unsupervised_auprc"]
    )
    result.to_csv(ROOT / "results" / "repeated_phyllosphere_transfer.csv", index=False)
    pd.concat(predictions, ignore_index=True).to_csv(
        ROOT / "analysis_data" / "repeated_phyllosphere_transfer_predictions.csv", index=False
    )
    agreement.to_csv(ROOT / "results" / "repeated_phyllosphere_label_agreement.csv")
    pd.DataFrame([{
        "experimentally_common_pairs": len(experimentally_common),
        "common_pairs_with_score_rows": len(common_pairs),
        "pairs_lost_to_score_coverage": len(experimentally_common - common_pairs),
        "label_agreements": int((labels["interaction_label_carlstrom"] == labels["interaction_label_schafer"]).sum()),
        "label_disagreements": int((labels["interaction_label_carlstrom"] != labels["interaction_label_schafer"]).sum()),
        "positive_in_both": int(((labels["interaction_label_carlstrom"] == 1) & (labels["interaction_label_schafer"] == 1)).sum()),
    }]).to_csv(ROOT / "results" / "repeated_phyllosphere_transfer_audit.csv", index=False)
    with open(ROOT / "results" / "repeated_phyllosphere_transfer_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "studies": STUDIES,
            "common_pairs": len(common_pairs),
            "experimentally_common_pairs": len(experimentally_common),
            "features": FEATURES,
            "validation": "train on one complete study and test on the other",
            "important_difference": (
                "Carlstrom measures removal effects; Schafer measures addition effects. "
                "The organisms overlap, but the biological endpoints differ."
            ),
        }, handle, indent=2)


if __name__ == "__main__":
    main()
