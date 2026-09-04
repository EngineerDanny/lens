#!/usr/bin/env python3
"""Compare network scores with direct pair summaries from abundance data."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    f1_score,
    precision_recall_curve,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from evaluate_within_system_calibration import FEATURES as SCORE_FEATURES
from evaluate_within_system_calibration import ROOT, SEED, UNSUPERVISED, load_system


SYSTEMS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
TAXON_FRACTIONS = [0.20, 0.40, 0.60]
REPEATS = 200
DIRECT_FEATURES = [
    "direct_prevalence_min",
    "direct_prevalence_max",
    "direct_joint_prevalence",
    "direct_presence_jaccard",
    "presence_phi_abs",
    "presence_log_odds_abs",
    "spearman_all_abs",
    "spearman_copresent_abs",
    "log_ratio_variance",
    "abundance_contrast_min",
    "abundance_contrast_max",
]
FEATURE_SETS = {
    "network_scores": SCORE_FEATURES,
    "direct_abundance": DIRECT_FEATURES,
    "combined": SCORE_FEATURES + DIRECT_FEATURES,
}


def safe_corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan
    return float(pd.Series(x).corr(pd.Series(y), method="spearman"))


def standardized_contrast(values: np.ndarray, group: np.ndarray) -> float:
    if group.sum() < 2 or (~group).sum() < 2:
        return np.nan
    pooled = np.sqrt((np.var(values[group], ddof=1) + np.var(values[~group], ddof=1)) / 2)
    if not np.isfinite(pooled) or pooled == 0:
        return 0.0
    return float((np.mean(values[group]) - np.mean(values[~group])) / pooled)


def build_direct_features(system: str) -> pd.DataFrame:
    abundance = pd.read_csv(ROOT / "cleaned_data" / f"{system}_abundance.csv")
    abundance = abundance.set_index(abundance.columns[0])
    x = abundance.apply(pd.to_numeric, errors="raise")
    taxa = list(x.columns)
    values = x.to_numpy(float)
    presence = values > 0
    n = len(x)
    positive = values[values > 0]
    pseudocount = float(np.min(positive) / 2) if positive.size else 0.5
    log_values = np.log(values + pseudocount)
    rows = []
    for first_index in range(len(taxa)):
        for second_index in range(first_index + 1, len(taxa)):
            first, second = taxa[first_index], taxa[second_index]
            p1, p2 = presence[:, first_index], presence[:, second_index]
            both = int(np.sum(p1 & p2))
            only1 = int(np.sum(p1 & ~p2))
            only2 = int(np.sum(~p1 & p2))
            neither = int(np.sum(~p1 & ~p2))
            union = both + only1 + only2
            prev1, prev2 = float(p1.mean()), float(p2.mean())
            phi = (np.corrcoef(p1.astype(float), p2.astype(float))[0, 1]
                   if np.unique(p1).size > 1 and np.unique(p2).size > 1 else 0.0)
            log_odds = np.log(((both + 0.5) * (neither + 0.5))
                              / ((only1 + 0.5) * (only2 + 0.5)))
            spearman_all = safe_corr(values[:, first_index], values[:, second_index])
            copresent = p1 & p2
            spearman_positive = safe_corr(
                values[copresent, first_index], values[copresent, second_index]
            )
            log_ratio_variance = float(np.var(
                log_values[:, first_index] - log_values[:, second_index], ddof=1
            ))
            first_contrast = standardized_contrast(log_values[:, first_index], p2)
            second_contrast = standardized_contrast(log_values[:, second_index], p1)
            contrasts = np.abs([first_contrast, second_contrast])
            rows.append({
                "taxon_1": min(str(first), str(second)),
                "taxon_2": max(str(first), str(second)),
                "direct_prevalence_min": min(prev1, prev2),
                "direct_prevalence_max": max(prev1, prev2),
                "direct_joint_prevalence": both / n,
                "direct_presence_jaccard": both / union if union else 0.0,
                "presence_phi_abs": abs(phi) if np.isfinite(phi) else 0.0,
                "presence_log_odds_abs": abs(log_odds),
                "spearman_all_abs": abs(spearman_all) if np.isfinite(spearman_all) else np.nan,
                "spearman_copresent_abs": abs(spearman_positive)
                if np.isfinite(spearman_positive) else np.nan,
                "log_ratio_variance": log_ratio_variance,
                "abundance_contrast_min": float(np.nanmin(contrasts))
                if np.isfinite(contrasts).any() else np.nan,
                "abundance_contrast_max": float(np.nanmax(contrasts))
                if np.isfinite(contrasts).any() else np.nan,
            })
    return pd.DataFrame(rows)


def load_analysis(system: str) -> pd.DataFrame:
    frame = load_system(system)
    direct = build_direct_features(system)
    frame = frame.merge(direct, on=["taxon_1", "taxon_2"], how="left", validate="many_to_one")
    if frame[DIRECT_FEATURES].isna().all(axis=0).any():
        missing = frame[DIRECT_FEATURES].columns[frame[DIRECT_FEATURES].isna().all()].tolist()
        raise ValueError(f"Entire direct feature columns missing for {system}: {missing}")
    return frame


def make_model() -> Pipeline:
    return Pipeline([
        ("impute", SimpleImputer(strategy="median", keep_empty_features=True)),
        ("scale", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l2", solver="liblinear", C=0.01, class_weight="balanced",
            max_iter=5000, random_state=SEED,
        )),
    ])


def score_model(calibration: pd.DataFrame, test: pd.DataFrame,
                features: list[str]) -> tuple[np.ndarray, np.ndarray, bool]:
    if calibration.interaction_label.nunique() < 2:
        probability = float(calibration.interaction_label.mean())
        return (np.repeat(probability, len(calibration)),
                np.repeat(probability, len(test)), False)
    model = make_model()
    model.fit(calibration[features], calibration.interaction_label)
    return (model.predict_proba(calibration[features])[:, 1],
            model.predict_proba(test[features])[:, 1], True)


def calibration_threshold(labels: pd.Series, scores: np.ndarray) -> float:
    """Maximize calibration F1; break ties toward a higher threshold."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=int)
    if labels.sum() == 0:
        return float(np.nextafter(scores.max(), np.inf))
    precision, recall, thresholds = precision_recall_curve(labels, scores)
    denominator = precision[:-1] + recall[:-1]
    values = np.divide(
        2 * precision[:-1] * recall[:-1], denominator,
        out=np.zeros_like(denominator), where=denominator > 0,
    )
    best = np.flatnonzero(np.isclose(values, values.max()))[-1]
    return float(thresholds[best])


def select_unsupervised(calibration: pd.DataFrame) -> tuple[str, str, bool]:
    """Choose a method using revealed labels only; use mean rank for one class."""
    if calibration.interaction_label.nunique() < 2:
        return "Mean rank", "mean_score_rank", False
    candidates = []
    for method, column in UNSUPERVISED.items():
        candidates.append((
            average_precision_score(calibration.interaction_label, calibration[column]),
            method,
            column,
        ))
    _, method, column = max(candidates)
    return method, column, True


def evaluate(frame: pd.DataFrame, calibration_pairs: set[str], test_pairs: set[str],
             system: str, design: str, fraction: float, repeat: int,
             selected_taxa: int) -> list[dict]:
    calibration = frame[frame.pair_id.isin(calibration_pairs)]
    test = frame[frame.pair_id.isin(test_pairs)]
    if calibration.empty or test.empty or test.interaction_label.nunique() < 2:
        return []
    y = test.interaction_label.to_numpy()
    rows = []
    selected_method, selected_column, selected_from_labels = select_unsupervised(calibration)
    selected_threshold = calibration_threshold(
        calibration.interaction_label, calibration[selected_column].to_numpy()
    )
    selected_test_scores = test[selected_column].to_numpy()
    rows.append({
        "analysis_set": system,
        "design": design,
        "taxon_fraction": fraction,
        "repeat": repeat,
        "selected_taxa": selected_taxa,
        "calibration_pairs": len(calibration_pairs),
        "calibration_positives": int(calibration.interaction_label.sum()),
        "test_pairs": len(test_pairs),
        "test_positives": int(y.sum()),
        "feature_set": "selected_unsupervised",
        "model_fitted": selected_from_labels,
        "selected_method": selected_method,
        "auprc": average_precision_score(y, selected_test_scores),
        "auroc": roc_auc_score(y, selected_test_scores),
        "f1": f1_score(y, selected_test_scores >= selected_threshold, zero_division=0),
        "threshold": selected_threshold,
    })
    for feature_set, features in FEATURE_SETS.items():
        calibration_score, score, fitted = score_model(calibration, test, features)
        threshold = calibration_threshold(calibration.interaction_label, calibration_score)
        rows.append({
            "analysis_set": system,
            "design": design,
            "taxon_fraction": fraction,
            "repeat": repeat,
            "selected_taxa": selected_taxa,
            "calibration_pairs": len(calibration_pairs),
            "calibration_positives": int(calibration.interaction_label.sum()),
            "test_pairs": len(test_pairs),
            "test_positives": int(y.sum()),
            "feature_set": feature_set,
            "model_fitted": fitted,
            "selected_method": "",
            "auprc": average_precision_score(y, score),
            "auroc": roc_auc_score(y, score),
            "f1": f1_score(y, score >= threshold, zero_division=0),
            "threshold": threshold,
        })
    return rows


def run_system(system: str) -> tuple[list[dict], list[dict]]:
    frame = load_analysis(system)
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
            panel_rows = evaluate(frame, panel_pairs, unseen_pairs, system, "taxon_panel",
                                  fraction, repeat, panel_size)
            if not panel_rows:
                failures.append({"analysis_set": system, "design": "taxon_panel",
                                 "taxon_fraction": fraction, "repeat": repeat,
                                 "reason": "test set lacks both classes"})
            rows.extend(panel_rows)
            budget = len(panel_pairs)
            if budget >= len(pairs):
                continue
            chosen = set(rng.choice(pairs.pair_id.to_numpy(), size=budget, replace=False))
            distributed_rows = evaluate(
                frame, chosen, set(pairs.pair_id) - chosen, system, "distributed_pairs",
                fraction, repeat, panel_size,
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
        ["analysis_set", "design", "taxon_fraction", "feature_set"], as_index=False
    ).agg(
        valid_repeats=("repeat", "size"),
        mean_calibration_pairs=("calibration_pairs", "mean"),
        model_fit_rate=("model_fitted", "mean"),
        mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"),
        mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"),
        mean_f1=("f1", "mean"),
        sd_f1=("f1", "std"),
    )
    reference = summary[summary.feature_set.eq("network_scores")][
        ["analysis_set", "design", "taxon_fraction", "mean_auprc"]
    ].rename(columns={"mean_auprc": "score_only_auprc"})
    summary = summary.merge(reference, on=["analysis_set", "design", "taxon_fraction"])
    summary["change_from_scores"] = summary.mean_auprc - summary.score_only_auprc
    results.to_csv(ROOT / "results" / "direct_pair_feature_repeats.csv", index=False)
    summary.to_csv(ROOT / "results" / "direct_pair_feature_summary.csv", index=False)
    pd.DataFrame(failures, columns=[
        "analysis_set", "design", "taxon_fraction", "repeat", "reason"
    ]).to_csv(ROOT / "results" / "direct_pair_feature_failures.csv", index=False)
    with open(ROOT / "results" / "direct_pair_feature_metadata.json", "w") as handle:
        json.dump({
            "seed": SEED,
            "repeats": REPEATS,
            "taxon_fractions": TAXON_FRACTIONS,
            "feature_sets": FEATURE_SETS,
            "unsupervised_selection": "highest calibration AUPRC among PLNNetwork, Poisson GLMNet, SparCC, and mean rank; mean rank fallback for one-class calibration sets",
            "threshold_selection": "threshold maximizing F1 on revealed calibration labels; ties resolved toward the higher threshold",
            "selection_uses_labels": False,
            "model": "fixed L2 logistic regression, C=0.01, balanced class weights",
            "pseudocount": "half the smallest positive abundance within each system",
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
