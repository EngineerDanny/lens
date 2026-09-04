#!/usr/bin/env python3
"""Evaluate sign-aware network inference results.

The sign label is defined only for experimentally tested pairs:

* tested neutral rows are labeled neutral;
* positive rows with a positive effect sign are labeled positive;
* positive rows with a negative effect sign are labeled negative;
* ambiguous rows and unsigned positives are excluded.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_fscore_support,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
CLEANED = ROOT / "cleaned_data"
FEATURES = ROOT / "analysis_data"
RESULTS = ROOT / "results"

METHODS = ["pln", "glmnet", "spieceasi", "spring", "sparcc"]
CLASSES = ["negative", "neutral", "positive"]


def normalize_sign(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"1", "1.0", "positive", "pos", "+", "+1"}:
        return "positive"
    if text in {"-1", "-1.0", "negative", "neg", "-"}:
        return "negative"
    if text in {"0", "0.0", "neutral", "none", "no_effect"}:
        return "neutral"
    return None


def normalize_status(value: object) -> str | None:
    if pd.isna(value):
        return None
    text = str(value).strip().lower()
    if text in {"neutral", "tested_neutral"}:
        return "neutral"
    if text == "positive":
        return "positive"
    if text == "ambiguous":
        return "ambiguous"
    return text


def derive_sign_label(row: pd.Series) -> str | None:
    if row["status_norm"] == "ambiguous":
        return None
    if row["status_norm"] == "neutral":
        return "neutral"
    if row["status_norm"] == "positive" and row["sign_norm"] in {
        "positive",
        "negative",
    }:
        return row["sign_norm"]
    return None


def load_sign_table() -> pd.DataFrame:
    parts: list[pd.DataFrame] = []
    for truth_path in sorted(CLEANED.glob("*_tested_pairs.csv")):
        dataset = truth_path.name.removesuffix("_tested_pairs.csv")
        feature_path = FEATURES / f"{dataset}_pair_features.csv"
        if not feature_path.exists():
            continue

        truth = pd.read_csv(truth_path)
        truth["status_norm"] = truth["tested_status"].map(normalize_status)
        truth["sign_norm"] = truth["effect_sign"].map(normalize_sign)
        truth["sign_label"] = truth.apply(derive_sign_label, axis=1)
        truth = truth[truth["sign_label"].notna()].copy()
        if truth.empty:
            continue

        features = pd.read_csv(feature_path)
        merged = truth.merge(
            features,
            on=["dataset", "taxon_1", "taxon_2"],
            how="inner",
            validate="one_to_one",
        )
        parts.append(merged)

    if not parts:
        raise RuntimeError("No usable signed tested pairs found.")
    return pd.concat(parts, ignore_index=True)


def metric_row(method: str, y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, object]:
    _, _, f1, _ = precision_recall_fscore_support(
        y_true,
        y_pred,
        labels=CLASSES,
        zero_division=0,
    )
    return {
        "method": method,
        "n": len(y_true),
        "pred_negative": int((y_pred == "negative").sum()),
        "pred_neutral": int((y_pred == "neutral").sum()),
        "pred_positive": int((y_pred == "positive").sum()),
        "accuracy": accuracy_score(y_true, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, labels=CLASSES, average="macro", zero_division=0),
        "negative_f1": f1[0],
        "neutral_f1": f1[1],
        "positive_f1": f1[2],
    }


def native_classification(sign_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    y_true = sign_table["sign_label"].to_numpy()
    for method in METHODS:
        score = sign_table[f"{method}_score"].fillna(0).to_numpy()
        selected = sign_table[f"{method}_selected"].fillna(0).astype(bool).to_numpy()
        y_pred = np.where(
            selected & (score > 0),
            "positive",
            np.where(selected & (score < 0), "negative", "neutral"),
        )
        rows.append(metric_row(method, y_true, y_pred))
    return pd.DataFrame(rows).sort_values("macro_f1", ascending=False)


def signed_score_ap(sign_table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    y_true = sign_table["sign_label"].to_numpy()
    baselines = {label: float((y_true == label).mean()) for label in CLASSES}
    for method in METHODS:
        score = sign_table[f"{method}_score"].fillna(0).to_numpy()
        scores = {
            "negative": -score,
            "neutral": -np.abs(score),
            "positive": score,
        }
        aps = {}
        lifts = {}
        for label in CLASSES:
            indicator = (y_true == label).astype(int)
            aps[label] = average_precision_score(indicator, scores[label])
            lifts[label] = aps[label] / baselines[label]
        rows.append(
            {
                "method": method,
                "n": len(y_true),
                "negative_ap": aps["negative"],
                "neutral_ap": aps["neutral"],
                "positive_ap": aps["positive"],
                "macro_ap": float(np.mean(list(aps.values()))),
                "negative_ap_lift": lifts["negative"],
                "neutral_ap_lift": lifts["neutral"],
                "positive_ap_lift": lifts["positive"],
                "macro_ap_lift": float(np.mean(list(lifts.values()))),
            }
        )
    return pd.DataFrame(rows).sort_values("macro_ap", ascending=False)


def supervised_multinomial(sign_table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_cols = []
    for method in METHODS:
        for suffix in ["score", "abs_score", "score_percentile", "selected"]:
            column = f"{method}_{suffix}"
            if column in sign_table.columns:
                feature_cols.append(column)

    for column in [
        "prevalence_1",
        "prevalence_2",
        "prevalence_min",
        "prevalence_max",
        "joint_prevalence",
        "n_joint_positive",
        "n_taxon_1_only",
        "n_taxon_2_only",
        "n_neither",
        "presence_jaccard",
        "n_samples",
        "n_taxa",
        "overall_zero_frequency",
    ]:
        if column in sign_table.columns:
            feature_cols.append(column)

    predictions = []
    for dataset in sorted(sign_table["dataset"].unique()):
        train = sign_table[sign_table["dataset"] != dataset].copy()
        test = sign_table[sign_table["dataset"] == dataset].copy()
        if set(train["sign_label"].unique()) != set(CLASSES):
            continue

        model = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median")),
                ("scaler", StandardScaler()),
                (
                    "classifier",
                    LogisticRegression(
                        max_iter=5000,
                        C=0.3,
                        class_weight="balanced",
                        solver="lbfgs",
                        random_state=1,
                    ),
                ),
            ]
        )
        model.fit(train[feature_cols], train["sign_label"])
        probabilities = model.predict_proba(test[feature_cols])
        labels = list(model.named_steps["classifier"].classes_)
        y_pred = np.array(labels)[np.argmax(probabilities, axis=1)]
        out = test[["dataset", "taxon_1", "taxon_2", "sign_label"]].copy()
        out["predicted_label"] = y_pred
        for index, label in enumerate(labels):
            out[f"probability_{label}"] = probabilities[:, index]
        predictions.append(out)

    pred = pd.concat(predictions, ignore_index=True)
    y_true = pred["sign_label"].to_numpy()
    y_pred = pred["predicted_label"].to_numpy()
    classification = pd.DataFrame(
        [metric_row("supervised_multinomial_calibration", y_true, y_pred)]
    )

    baselines = {label: float((sign_table["sign_label"].to_numpy() == label).mean()) for label in CLASSES}
    aps = {
        label: average_precision_score(
            (y_true == label).astype(int),
            pred[f"probability_{label}"],
        )
        for label in CLASSES
    }
    ap = pd.DataFrame(
        [
            {
                "method": "supervised_multinomial_calibration",
                "n": len(y_true),
                "negative_ap": aps["negative"],
                "neutral_ap": aps["neutral"],
                "positive_ap": aps["positive"],
                "macro_ap": float(np.mean(list(aps.values()))),
                "negative_ap_lift": aps["negative"] / baselines["negative"],
                "neutral_ap_lift": aps["neutral"] / baselines["neutral"],
                "positive_ap_lift": aps["positive"] / baselines["positive"],
                "macro_ap_lift": float(np.mean([aps[label] / baselines[label] for label in CLASSES])),
            }
        ]
    )

    confusion = pd.DataFrame(
        confusion_matrix(y_true, y_pred, labels=CLASSES),
        index=[f"true_{label}" for label in CLASSES],
        columns=[f"pred_{label}" for label in CLASSES],
    )
    return classification, ap, confusion


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    sign_table = load_sign_table()

    label_counts = (
        pd.crosstab(sign_table["dataset"], sign_table["sign_label"])
        .reindex(columns=CLASSES, fill_value=0)
        .reset_index()
    )
    label_counts.to_csv(RESULTS / "sign_label_counts_by_dataset.csv", index=False)
    native_classification(sign_table).to_csv(
        RESULTS / "sign_native_classification.csv",
        index=False,
    )
    signed_score_ap(sign_table).to_csv(
        RESULTS / "sign_score_one_vs_rest_ap.csv",
        index=False,
    )
    classification, ap, confusion = supervised_multinomial(sign_table)
    classification.to_csv(
        RESULTS / "sign_supervised_multinomial_classification.csv",
        index=False,
    )
    ap.to_csv(RESULTS / "sign_supervised_multinomial_ap.csv", index=False)
    confusion.to_csv(RESULTS / "sign_supervised_multinomial_confusion.csv")

    print("Wrote sign analysis results to results/")


if __name__ == "__main__":
    main()
