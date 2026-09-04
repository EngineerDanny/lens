#!/usr/bin/env python3
"""Leakage-safe signed LENS analysis for the Carlstrom system.

The analysis uses canonical undirected taxon pairs. Experimental effects are
negative, neutral, or positive. Five stratified outer folds produce one
held-out prediction for every scored pair. A class-weighted multinomial ridge
model is tuned by three-fold stratified validation within each training fold.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, label_binarize


ROOT = Path(__file__).resolve().parents[1]
SYSTEM = "carlstrom_phyllosphere_2019"
SEED = 20260821
CLASSES = np.array(["negative", "neutral", "positive"])
C_VALUES = (0.001, 0.01, 0.1, 1.0, 10.0)
FEATURES = [
    "pln_score",
    "glmnet_score",
    "sparcc_score",
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
    "prevalence_min",
    "prevalence_max",
    "joint_prevalence",
    "presence_jaccard",
]


def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    first = frame["taxon_1"].astype(str)
    second = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(first, second)
    frame["taxon_2"] = np.maximum(first, second)
    frame["pair_id"] = frame["taxon_1"] + "||" + frame["taxon_2"]
    return frame


def sign_label(row: pd.Series) -> str | None:
    status = str(row["tested_status"]).strip().lower()
    sign = str(row["effect_sign"]).strip().lower()
    if status in {"neutral", "tested_neutral"}:
        return "neutral"
    if status == "positive" and sign in {"negative", "-1", "-1.0"}:
        return "negative"
    if status == "positive" and sign in {"positive", "1", "1.0"}:
        return "positive"
    return None


def load_data() -> tuple[pd.DataFrame, dict[str, object]]:
    truth = canonicalize(pd.read_csv(
        ROOT / "cleaned_data" / f"{SYSTEM}_tested_pairs.csv"
    ))
    features = canonicalize(pd.read_csv(
        ROOT / "analysis_data" / f"{SYSTEM}_pair_features.csv"
    ))
    onenet = canonicalize(pd.read_csv(
        ROOT / "analysis_data" / f"{SYSTEM}_onenet_pair_scores.csv"
    ))
    if truth["pair_id"].duplicated().any() or features["pair_id"].duplicated().any():
        raise ValueError("Canonical pair identifiers must be unique.")

    truth["sign_label"] = truth.apply(sign_label, axis=1)
    merged = truth.merge(
        features.drop(columns=["dataset", "taxon_1", "taxon_2"]),
        on="pair_id",
        how="inner",
        validate="one_to_one",
    )
    merged = merged.merge(
        onenet[["pair_id", "onenet_score"]],
        on="pair_id",
        how="left",
        validate="one_to_one",
    )
    merged = merged[merged["sign_label"].notna()].copy().reset_index(drop=True)
    if merged["onenet_score"].isna().any():
        raise ValueError("OneNet scores are missing after canonical alignment.")
    if merged[FEATURES].isna().all(axis=0).any():
        raise ValueError("At least one signed LENS feature is entirely missing.")

    feature_taxa = set(features["taxon_1"]) | set(features["taxon_2"])
    unmatched = truth[~truth["pair_id"].isin(features["pair_id"])]
    audit = {
        "truth_pairs": int(len(truth)),
        "feature_pairs": int(len(features)),
        "canonical_matches": int(len(merged)),
        "unmatched_pairs": int(len(unmatched)),
        "unmatched_with_both_taxa_present": int((
            unmatched["taxon_1"].isin(feature_taxa)
            & unmatched["taxon_2"].isin(feature_taxa)
        ).sum()),
        "truth_only_taxa": sorted(
            (set(truth["taxon_1"]) | set(truth["taxon_2"])) - feature_taxa
        ),
        "class_counts": {
            label: int((merged["sign_label"] == label).sum()) for label in CLASSES
        },
    }
    return merged, audit


def make_model(c_value: float) -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2",
            solver="lbfgs",
            C=c_value,
            class_weight="balanced",
            max_iter=10_000,
            random_state=SEED,
        )),
    ])


def macro_ap(y: np.ndarray, probability: np.ndarray) -> tuple[float, dict[str, float]]:
    binary = label_binarize(y, classes=CLASSES)
    values = {
        label: float(average_precision_score(binary[:, index], probability[:, index]))
        for index, label in enumerate(CLASSES)
    }
    return float(np.mean(list(values.values()))), values


def tune_c(train: pd.DataFrame, seed: int) -> tuple[float, pd.DataFrame]:
    splitter = StratifiedKFold(n_splits=3, shuffle=True, random_state=seed)
    rows: list[dict[str, float]] = []
    y = train["sign_label"].to_numpy()
    for c_value in C_VALUES:
        probability = np.full((len(train), len(CLASSES)), np.nan)
        for fit_index, validation_index in splitter.split(train, y):
            model = make_model(c_value)
            model.fit(train.iloc[fit_index][FEATURES], y[fit_index])
            fold_probability = model.predict_proba(train.iloc[validation_index][FEATURES])
            order = [list(model.named_steps["model"].classes_).index(x) for x in CLASSES]
            probability[validation_index] = fold_probability[:, order]
        score, per_class = macro_ap(y, probability)
        rows.append({"C": c_value, "macro_ap": score, **{
            f"{label}_ap": per_class[label] for label in CLASSES
        }})
    tuning = pd.DataFrame(rows)
    best = tuning.sort_values(["macro_ap", "C"], ascending=[False, True]).iloc[0]
    return float(best["C"]), tuning


def signed_unsupervised_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    y = frame["sign_label"].to_numpy()
    for method in ("pln", "glmnet", "sparcc"):
        score = frame[f"{method}_score"].fillna(0).to_numpy(float)
        evidence = np.column_stack((-score, -np.abs(score), score))
        per_class = {}
        for index, label in enumerate(CLASSES):
            per_class[label] = average_precision_score(y == label, evidence[:, index])
        rows.append({
            "method": method,
            **{f"{label}_ap": per_class[label] for label in CLASSES},
            "macro_ap": float(np.mean(list(per_class.values()))),
        })
    return pd.DataFrame(rows)


def main() -> None:
    frame, audit = load_data()
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    all_predictions = []
    tuning_rows = []
    fold_rows = []

    for fold, (train_index, test_index) in enumerate(
        outer.split(frame, frame["sign_label"]), start=1
    ):
        train = frame.iloc[train_index].copy()
        test = frame.iloc[test_index].copy()
        selected_c, tuning = tune_c(train, SEED + fold)
        tuning["outer_fold"] = fold
        tuning["selected"] = tuning["C"] == selected_c
        tuning_rows.append(tuning)

        model = make_model(selected_c)
        model.fit(train[FEATURES], train["sign_label"])
        raw_probability = model.predict_proba(test[FEATURES])
        model_classes = list(model.named_steps["model"].classes_)
        order = [model_classes.index(label) for label in CLASSES]
        probability = raw_probability[:, order]
        predicted = CLASSES[np.argmax(probability, axis=1)]
        interaction_probability = 1.0 - probability[:, 1]

        # The display network uses a training-only estimate of interaction
        # prevalence to determine how many held-out pairs receive an edge.
        training_prevalence = float((train["sign_label"] != "neutral").mean())
        display_edges = max(1, int(round(training_prevalence * len(test))))
        selected = np.zeros(len(test), dtype=bool)
        selected[np.argsort(-interaction_probability)[:display_edges]] = True
        onenet_selected = np.zeros(len(test), dtype=bool)
        onenet_selected[
            np.argsort(-test["onenet_score"].to_numpy(float))[:display_edges]
        ] = True
        display_sign = np.where(
            probability[:, 0] >= probability[:, 2], "negative", "positive"
        )

        out = test[[
            "pair_id", "taxon_1", "taxon_2", "sign_label",
            "prevalence_1", "prevalence_2",
            "onenet_score",
        ]].copy()
        out["fold"] = fold
        out["selected_C"] = selected_c
        out["predicted_label"] = predicted
        out["probability_negative"] = probability[:, 0]
        out["probability_neutral"] = probability[:, 1]
        out["probability_positive"] = probability[:, 2]
        out["interaction_probability"] = interaction_probability
        out["display_selected"] = selected
        out["display_sign"] = display_sign
        out["onenet_display_selected"] = onenet_selected
        all_predictions.append(out)

        fold_macro, fold_ap = macro_ap(test["sign_label"].to_numpy(), probability)
        fold_rows.append({
            "fold": fold,
            "train_pairs": len(train),
            "test_pairs": len(test),
            "train_negative": int((train["sign_label"] == "negative").sum()),
            "train_neutral": int((train["sign_label"] == "neutral").sum()),
            "train_positive": int((train["sign_label"] == "positive").sum()),
            "selected_C": selected_c,
            "negative_ap": fold_ap["negative"],
            "neutral_ap": fold_ap["neutral"],
            "positive_ap": fold_ap["positive"],
            "macro_ap": fold_macro,
        })

    predictions = pd.concat(all_predictions, ignore_index=True)
    probability = predictions[[f"probability_{x}" for x in CLASSES]].to_numpy()
    y = predictions["sign_label"].to_numpy()
    predicted = predictions["predicted_label"].to_numpy()
    overall_macro, overall_ap = macro_ap(y, probability)
    overall = pd.DataFrame([{
        "method": "signed_LENS",
        "pairs": len(predictions),
        "negative_pairs": int((y == "negative").sum()),
        "neutral_pairs": int((y == "neutral").sum()),
        "positive_pairs": int((y == "positive").sum()),
        "negative_ap": overall_ap["negative"],
        "neutral_ap": overall_ap["neutral"],
        "positive_ap": overall_ap["positive"],
        "macro_ap": overall_macro,
        "balanced_accuracy": balanced_accuracy_score(y, predicted),
        "macro_f1": f1_score(y, predicted, labels=CLASSES, average="macro"),
    }])
    baselines = signed_unsupervised_metrics(frame)
    comparison = pd.concat([
        overall[["method", "negative_ap", "neutral_ap", "positive_ap", "macro_ap"]],
        baselines,
    ], ignore_index=True)
    confusion = pd.DataFrame(
        confusion_matrix(y, predicted, labels=CLASSES),
        index=[f"true_{x}" for x in CLASSES],
        columns=[f"pred_{x}" for x in CLASSES],
    )
    interaction = (predictions["sign_label"] != "neutral").astype(int)
    lens_selected = predictions["display_selected"].astype(bool)
    onenet_selected = predictions["onenet_display_selected"].astype(bool)
    edge_summary = pd.DataFrame([
        {
            "method": "LENS",
            "selected_edges": int(lens_selected.sum()),
            "selected_positive_effects": int(((predictions["sign_label"] == "positive") & lens_selected).sum()),
            "selected_negative_effects": int(((predictions["sign_label"] == "negative") & lens_selected).sum()),
            "selected_neutrals": int(((predictions["sign_label"] == "neutral") & lens_selected).sum()),
            "precision": float(interaction[lens_selected].mean()),
            "recall": float(interaction[lens_selected].sum() / interaction.sum()),
            "interaction_auprc": float(average_precision_score(interaction, predictions["interaction_probability"])),
        },
        {
            "method": "OneNet",
            "selected_edges": int(onenet_selected.sum()),
            "selected_positive_effects": int(((predictions["sign_label"] == "positive") & onenet_selected).sum()),
            "selected_negative_effects": int(((predictions["sign_label"] == "negative") & onenet_selected).sum()),
            "selected_neutrals": int(((predictions["sign_label"] == "neutral") & onenet_selected).sum()),
            "precision": float(interaction[onenet_selected].mean()),
            "recall": float(interaction[onenet_selected].sum() / interaction.sum()),
            "interaction_auprc": float(average_precision_score(interaction, predictions["onenet_score"])),
        },
    ])

    results = ROOT / "results"
    analysis = ROOT / "analysis_data"
    results.mkdir(exist_ok=True)
    predictions.to_csv(analysis / "carlstrom_signed_lens_oof_predictions.csv", index=False)
    pd.concat(tuning_rows, ignore_index=True).to_csv(
        results / "carlstrom_signed_lens_inner_tuning.csv", index=False
    )
    pd.DataFrame(fold_rows).to_csv(
        results / "carlstrom_signed_lens_fold_metrics.csv", index=False
    )
    overall.to_csv(results / "carlstrom_signed_lens_metrics.csv", index=False)
    comparison.to_csv(results / "carlstrom_signed_lens_comparison.csv", index=False)
    confusion.to_csv(results / "carlstrom_signed_lens_confusion.csv")
    edge_summary.to_csv(results / "carlstrom_signed_network_edge_summary.csv", index=False)
    with open(results / "carlstrom_signed_lens_audit.json", "w") as handle:
        json.dump({
            **audit,
            "seed": SEED,
            "outer_validation": "five-fold stratified cross-validation",
            "inner_validation": "three-fold stratified cross-validation",
            "model": "class-weighted multinomial ridge logistic regression",
            "features": FEATURES,
            "display_edge_budget": "training interaction prevalence applied to each held-out fold",
        }, handle, indent=2)

    print("Alignment audit")
    print(json.dumps(audit, indent=2))
    print("\nHeld-out performance")
    print(overall.to_string(index=False))
    print("\nComparison")
    print(comparison.sort_values("macro_ap", ascending=False).to_string(index=False))
    print("\nConfusion matrix")
    print(confusion.to_string())
    print("\nNetwork edge comparison")
    print(edge_summary.to_string(index=False))


if __name__ == "__main__":
    main()
