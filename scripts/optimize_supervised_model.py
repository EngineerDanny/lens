#!/usr/bin/env python3
"""Nested optimization of supervised pair ranking for AUPRC."""

from __future__ import annotations

import json
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.exceptions import ConvergenceWarning
from scipy.optimize import minimize
from scipy.special import expit

from evaluate_direct_pair_features import FEATURE_SETS
from evaluate_pair_five_fold_cv import (
    BUDGET_FRACTIONS,
    REFERENCE_METHODS,
    SYSTEMS,
    load_system,
)
from evaluate_within_system_calibration import ROOT, SEED


FEATURES = FEATURE_SETS["combined"]
FIXED_FALLBACK_FEATURES = [
    "pln_score_percentile",
    "glmnet_score_percentile",
    "sparcc_score_percentile",
]
MAX_RANKING_PAIRS = 30_000


@dataclass(frozen=True)
class Candidate:
    family: str
    setting: str


CANDIDATES = [
    Candidate("logistic", "l2_C=0.01"),
    Candidate("logistic", "l2_C=0.1"),
    Candidate("logistic", "l2_C=1"),
    Candidate("logistic", "l1_C=0.01"),
    Candidate("logistic", "l1_C=0.1"),
    Candidate("logistic", "l1_C=1"),
    Candidate("ranking", "C=0.01"),
    Candidate("ranking", "C=0.1"),
    Candidate("ranking", "C=1"),
    Candidate("boosting", "depth=1_l2=1"),
    Candidate("boosting", "depth=2_l2=1"),
]

warnings.filterwarnings("ignore", category=FutureWarning, module="sklearn")
warnings.filterwarnings("ignore", category=ConvergenceWarning)
warnings.filterwarnings(
    "ignore", category=UserWarning, module="sklearn.linear_model._logistic"
)


def parse_value(setting: str, key: str) -> float:
    token = next(x for x in setting.split("_") if x.startswith(f"{key}="))
    return float(token.split("=", 1)[1])


def preprocess(train: pd.DataFrame, test: pd.DataFrame):
    imputer = SimpleImputer(strategy="median", keep_empty_features=True)
    scaler = StandardScaler()
    x_train = scaler.fit_transform(imputer.fit_transform(train[FEATURES]))
    x_test = scaler.transform(imputer.transform(test[FEATURES]))
    return x_train, x_test


def fit_score(candidate: Candidate, train: pd.DataFrame, test: pd.DataFrame,
              seed: int) -> np.ndarray:
    y = train["interaction_label"].to_numpy(int)
    class_counts = np.bincount(y, minlength=2)
    if class_counts.min() < 3:
        return test[FIXED_FALLBACK_FEATURES].mean(axis=1).to_numpy(float)
    x_train, x_test = preprocess(train, test)

    if candidate.family == "centered_ridge":
        # On the standardized scale, these coefficients reproduce the fixed
        # mean of the three raw score percentiles up to an intercept shift.
        imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        scaler = StandardScaler()
        x_train = scaler.fit_transform(imputer.fit_transform(train[FEATURES]))
        x_test = scaler.transform(imputer.transform(test[FEATURES]))
        reference = np.zeros(x_train.shape[1])
        for feature in FIXED_FALLBACK_FEATURES:
            feature_index = FEATURES.index(feature)
            reference[feature_index] = scaler.scale_[feature_index] / 3.0
        event_count = y.sum()
        nonevent_count = len(y) - event_count
        sample_weight = np.where(
            y == 1,
            len(y) / (2.0 * event_count),
            len(y) / (2.0 * nonevent_count),
        )
        penalty = parse_value(candidate.setting, "lambda")

        def objective(parameters: np.ndarray):
            intercept = parameters[0]
            coefficients = parameters[1:]
            linear_predictor = intercept + x_train @ coefficients
            loss = np.average(
                np.logaddexp(0.0, linear_predictor) - y * linear_predictor,
                weights=sample_weight,
            )
            displacement = coefficients - reference
            value = loss + 0.5 * penalty * np.dot(displacement, displacement)
            probability = expit(linear_predictor)
            residual = sample_weight * (probability - y)
            residual /= sample_weight.sum()
            gradient = np.r_[
                residual.sum(),
                x_train.T @ residual + penalty * displacement,
            ]
            return value, gradient

        initial = np.r_[0.0, reference]
        fitted = minimize(
            objective,
            initial,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 2_000, "ftol": 1e-12},
        )
        if not fitted.success:
            raise RuntimeError(f"Centered ridge failed: {fitted.message}")
        return expit(fitted.x[0] + x_test @ fitted.x[1:])

    if candidate.family == "logistic":
        c_value = parse_value(candidate.setting, "C")
        if candidate.setting.startswith("l1_"):
            model = LogisticRegression(
                penalty="l1", solver="liblinear", C=c_value,
                class_weight="balanced", max_iter=10_000,
                random_state=seed, tol=1e-6,
            )
        else:
            model = LogisticRegression(
                penalty="l2", solver="liblinear", C=c_value,
                class_weight="balanced", max_iter=5_000,
                random_state=seed,
            )
        model.fit(x_train, y)
        return model.predict_proba(x_test)[:, 1]

    if candidate.family == "ranking":
        positive = np.flatnonzero(y == 1)
        negative = np.flatnonzero(y == 0)
        total = len(positive) * len(negative)
        rng = np.random.default_rng(seed)
        if total <= MAX_RANKING_PAIRS:
            p = np.repeat(positive, len(negative))
            n = np.tile(negative, len(positive))
        else:
            p = rng.choice(positive, MAX_RANKING_PAIRS, replace=True)
            n = rng.choice(negative, MAX_RANKING_PAIRS, replace=True)
        differences = x_train[p] - x_train[n]
        # Include reversed comparisons to prevent an intercept from encoding order.
        rank_x = np.vstack([differences, -differences])
        rank_y = np.r_[np.ones(len(differences)), np.zeros(len(differences))]
        model = LogisticRegression(
            penalty="l2", solver="liblinear",
            C=parse_value(candidate.setting, "C"), fit_intercept=False,
            max_iter=5_000, random_state=seed,
        )
        model.fit(rank_x, rank_y)
        return x_test @ model.coef_[0]

    depth = int(parse_value(candidate.setting, "depth"))
    model = HistGradientBoostingClassifier(
        learning_rate=0.05, max_iter=100, max_depth=depth,
        min_samples_leaf=max(5, min(20, len(train) // 10)),
        l2_regularization=parse_value(candidate.setting, "l2"),
        random_state=seed,
    )
    class_weight = np.where(y == 1, len(y) / (2 * y.sum()),
                            len(y) / (2 * (len(y) - y.sum())))
    model.fit(x_train, y, sample_weight=class_weight)
    return model.predict_proba(x_test)[:, 1]


def grouped_folds(frame: pd.DataFrame, seed: int, maximum: int = 3):
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    class_counts = groups["group_label"].value_counts()
    folds = min(maximum, int(class_counts.min())) if len(class_counts) == 2 else 0
    if folds < 2:
        return []
    splitter = StratifiedKFold(n_splits=folds, shuffle=True, random_state=seed)
    output = []
    for train_idx, valid_idx in splitter.split(groups["pair_id"], groups["group_label"]):
        train_ids = set(groups.iloc[train_idx]["pair_id"])
        valid_ids = set(groups.iloc[valid_idx]["pair_id"])
        output.append((frame[frame.pair_id.isin(train_ids)],
                       frame[frame.pair_id.isin(valid_ids)]))
    return output


def select_candidate(training: pd.DataFrame, seed: int):
    folds = grouped_folds(training, seed)
    if not folds:
        return Candidate("logistic", "l2_C=0.01"), pd.DataFrame()
    rows = []
    for candidate_index, candidate in enumerate(CANDIDATES):
        values = []
        for fold_index, (inner_train, inner_valid) in enumerate(folds):
            y_valid = inner_valid.interaction_label.to_numpy(int)
            if np.unique(y_valid).size < 2:
                continue
            score = fit_score(candidate, inner_train, inner_valid,
                              seed + 1000 * candidate_index + fold_index)
            values.append(average_precision_score(y_valid, score))
        rows.append({
            "family": candidate.family,
            "setting": candidate.setting,
            "inner_folds": len(values),
            "mean_inner_auprc": float(np.mean(values)) if values else np.nan,
        })
    tuning = pd.DataFrame(rows)
    valid = tuning.dropna(subset=["mean_inner_auprc"])
    if valid.empty:
        return Candidate("logistic", "l2_C=0.01"), tuning
    best = valid.sort_values(
        ["mean_inner_auprc", "family", "setting"], ascending=[False, True, True]
    ).iloc[0]
    return Candidate(best.family, best.setting), tuning


def run_system(system: str):
    frame = load_system(system)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    result_rows, tuning_rows = [], []
    _, reference_column = REFERENCE_METHODS[system]

    for fold, (pool_index, test_index) in enumerate(
        outer.split(groups["pair_id"], groups["group_label"]), start=1
    ):
        pool_ids = groups.iloc[pool_index]["pair_id"].to_numpy()
        test_ids = set(groups.iloc[test_index]["pair_id"])
        test = frame[frame.pair_id.isin(test_ids)]
        y_test = test.interaction_label.to_numpy(int)
        rng = np.random.default_rng(SEED + 100_000 * SYSTEMS.index(system) + fold)
        ordered_pool = rng.permutation(pool_ids)
        for fraction in BUDGET_FRACTIONS:
            budget = min(int(np.ceil(fraction * len(groups))), len(ordered_pool))
            training_ids = set(ordered_pool[:budget])
            training = frame[frame.pair_id.isin(training_ids)]
            choice, tuning = select_candidate(
                training, SEED + 1_000_000 * SYSTEMS.index(system) + 10_000 * fold
                + int(100 * fraction)
            )
            score = fit_score(choice, training, test, SEED + fold)
            current = fit_score(Candidate("logistic", "l2_C=0.01"), training, test,
                                SEED + fold)
            result_rows.append({
                "analysis_set": system, "fold": fold,
                "budget_fraction": fraction, "budget_pairs": budget,
                "test_pairs": len(test_ids), "test_positives": int(y_test.sum()),
                "selected_family": choice.family, "selected_setting": choice.setting,
                "optimized_auprc": average_precision_score(y_test, score),
                "current_auprc": average_precision_score(y_test, current),
                "onenet_auprc": average_precision_score(y_test, test[reference_column]),
                "optimized_auroc": roc_auc_score(y_test, score),
                "current_auroc": roc_auc_score(y_test, current),
            })
            if not tuning.empty:
                tuning.insert(0, "analysis_set", system)
                tuning.insert(1, "fold", fold)
                tuning.insert(2, "budget_fraction", fraction)
                tuning_rows.extend(tuning.to_dict("records"))
    return result_rows, tuning_rows


def main():
    results, tuning = [], []
    for system in SYSTEMS:
        print(f"Optimizing {system}", flush=True)
        system_results, system_tuning = run_system(system)
        results.extend(system_results)
        tuning.extend(system_tuning)
    result = pd.DataFrame(results)
    summary = result.groupby(["analysis_set", "budget_fraction"], as_index=False).agg(
        budget_pairs=("budget_pairs", "max"), folds=("fold", "size"),
        current_auprc=("current_auprc", "mean"),
        optimized_auprc=("optimized_auprc", "mean"),
        onenet_auprc=("onenet_auprc", "mean"),
        optimized_sd=("optimized_auprc", "std"),
        current_sd=("current_auprc", "std"),
        optimized_auroc=("optimized_auroc", "mean"),
        current_auroc=("current_auroc", "mean"),
    )
    summary["absolute_improvement"] = (
        summary.optimized_auprc - summary.current_auprc
    )
    result.to_csv(ROOT / "results" / "supervised_model_optimization.csv", index=False)
    summary.to_csv(ROOT / "results" / "supervised_model_optimization_summary.csv", index=False)
    pd.DataFrame(tuning).to_csv(
        ROOT / "results" / "supervised_model_inner_tuning.csv", index=False
    )
    with open(ROOT / "results" / "supervised_model_optimization_metadata.json", "w") as handle:
        json.dump({
            "outer_folds": 5, "inner_folds_maximum": 3,
            "selection_metric": "mean inner-fold AUPRC",
            "pair_grouping": True, "features": FEATURES,
            "candidate_models": [candidate.__dict__ for candidate in CANDIDATES],
            "maximum_sampled_positive_negative_comparisons": MAX_RANKING_PAIRS,
            "test_labels_used_for_selection": False,
        }, handle, indent=2)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
