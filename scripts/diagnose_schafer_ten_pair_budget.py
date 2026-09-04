#!/usr/bin/env python3
"""Diagnose sparse logistic behavior at the 10-pair Schafer budget."""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

from evaluate_optimized_sparse_logistic import select_regularization
from evaluate_pair_five_fold_cv import REFERENCE_METHODS, SYSTEMS, load_system
from evaluate_within_system_calibration import ROOT, SEED
from optimize_supervised_model import FEATURES, parse_value


SYSTEM = "schafer_phyllosphere_2022"
REPEATS = 20
BUDGET = 10
SIMPLE_FEATURE_SETS = {
    "mean_three_scores": [
        "pln_score_percentile",
        "glmnet_score_percentile",
        "sparcc_score_percentile",
    ],
    "three_scores": [
        "pln_score_percentile",
        "glmnet_score_percentile",
        "sparcc_score_percentile",
    ],
    "three_scores_plus_joint_prevalence": [
        "pln_score_percentile",
        "glmnet_score_percentile",
        "sparcc_score_percentile",
        "direct_joint_prevalence",
    ],
}


def simple_score(
    name: str, features: list[str], training: pd.DataFrame,
    test: pd.DataFrame, seed: int,
) -> np.ndarray:
    y = training["interaction_label"].to_numpy(int)
    if np.unique(y).size < 2:
        return np.repeat(float(y.mean()), len(test))
    if name == "mean_three_scores":
        x_train = training[features].mean(axis=1).to_numpy().reshape(-1, 1)
        x_test = test[features].mean(axis=1).to_numpy().reshape(-1, 1)
    else:
        imputer = SimpleImputer(strategy="median", keep_empty_features=True)
        scaler = StandardScaler()
        x_train = scaler.fit_transform(imputer.fit_transform(training[features]))
        x_test = scaler.transform(imputer.transform(test[features]))
    model = LogisticRegression(
        penalty="l2", solver="liblinear", C=0.01,
        class_weight="balanced", max_iter=5_000, random_state=seed,
    )
    model.fit(x_train, y)
    return model.predict_proba(x_test)[:, 1]


def main() -> None:
    frame = load_system(SYSTEM)
    groups = frame.groupby("pair_id", as_index=False).agg(
        group_label=("interaction_label", "max")
    )
    _, onenet_column = REFERENCE_METHODS[SYSTEM]
    system_index = SYSTEMS.index(SYSTEM)
    rows = []

    for repeat in range(1, REPEATS + 1):
        outer_seed = SEED + 10_000_000 * repeat
        splitter = StratifiedKFold(n_splits=5, shuffle=True, random_state=outer_seed)
        for fold, (pool_index, test_index) in enumerate(
            splitter.split(groups["pair_id"], groups["group_label"]), start=1
        ):
            pool_ids = groups.iloc[pool_index]["pair_id"].to_numpy()
            test_ids = set(groups.iloc[test_index]["pair_id"])
            test = frame[frame["pair_id"].isin(test_ids)].copy()
            y_test = test["interaction_label"].to_numpy(int)

            order_seed = outer_seed + 100_000 * system_index + fold
            ordered_pool = np.random.default_rng(order_seed).permutation(pool_ids)
            training_ids = set(ordered_pool[:BUDGET])
            training = frame[frame["pair_id"].isin(training_ids)].copy()
            y_train = training["interaction_label"].to_numpy(int)
            candidate, tuning = select_regularization(training, order_seed)

            single_class = np.unique(y_train).size < 2
            if single_class:
                prediction = np.repeat(float(y_train.mean()), len(test))
                intercept = np.nan
                nonzero = 0
                coefficient_l1 = 0.0
            else:
                imputer = SimpleImputer(strategy="median", keep_empty_features=True)
                scaler = StandardScaler()
                x_train = scaler.fit_transform(imputer.fit_transform(training[FEATURES]))
                x_test = scaler.transform(imputer.transform(test[FEATURES]))
                model = LogisticRegression(
                    penalty="l1",
                    solver="liblinear",
                    C=parse_value(candidate.setting, "C"),
                    class_weight="balanced",
                    max_iter=10_000,
                    random_state=order_seed,
                    tol=1e-6,
                )
                model.fit(x_train, y_train)
                prediction = model.predict_proba(x_test)[:, 1]
                intercept = float(model.intercept_[0])
                nonzero = int(np.count_nonzero(model.coef_))
                coefficient_l1 = float(np.abs(model.coef_).sum())

            finite_tuning = [
                row["mean_inner_auprc"] for row in tuning
                if np.isfinite(row["mean_inner_auprc"])
            ]
            simple_predictions = {
                name: simple_score(name, features, training, test, order_seed)
                for name, features in SIMPLE_FEATURE_SETS.items()
            }
            fixed_mean_score = test[SIMPLE_FEATURE_SETS["mean_three_scores"]].mean(
                axis=1
            ).to_numpy(float)
            adaptive_prediction = (
                fixed_mean_score if single_class
                else simple_predictions["mean_three_scores"]
            )
            rows.append({
                "repeat": repeat,
                "fold": fold,
                "training_positives": int(y_train.sum()),
                "training_neutrals": int((y_train == 0).sum()),
                "single_class_training": single_class,
                "inner_tuning_available": bool(finite_tuning),
                "selected_setting": candidate.setting,
                "nonzero_coefficients": nonzero,
                "coefficient_l1_norm": coefficient_l1,
                "intercept": intercept,
                "prediction_min": float(prediction.min()),
                "prediction_max": float(prediction.max()),
                "prediction_sd": float(np.std(prediction)),
                "unique_predictions": int(np.unique(np.round(prediction, 12)).size),
                "test_pairs": len(test),
                "test_positives": int(y_test.sum()),
                "positive_frequency": float(y_test.mean()),
                "supervised_auprc": average_precision_score(y_test, prediction),
                "onenet_auprc": average_precision_score(y_test, test[onenet_column]),
                **{
                    f"{name}_auprc": average_precision_score(y_test, values)
                    for name, values in simple_predictions.items()
                },
                "fixed_mean_score_auprc": average_precision_score(
                    y_test, fixed_mean_score
                ),
                "adaptive_simple_auprc": average_precision_score(
                    y_test, adaptive_prediction
                ),
            })

    diagnostics = pd.DataFrame(rows)
    diagnostics.to_csv(
        ROOT / "results" / "schafer_ten_pair_diagnostics.csv", index=False
    )

    diagnostics["fit_state"] = np.select(
        [
            diagnostics.single_class_training,
            diagnostics.nonzero_coefficients.eq(0),
        ],
        ["no positive training example", "all coefficients zero"],
        default="nonzero fitted model",
    )
    state_summary = diagnostics.groupby("fit_state", as_index=False).agg(
        fits=("fold", "size"),
        mean_training_positives=("training_positives", "mean"),
        mean_prediction_sd=("prediction_sd", "mean"),
        mean_supervised_auprc=("supervised_auprc", "mean"),
        mean_onenet_auprc=("onenet_auprc", "mean"),
    )
    state_summary["fit_fraction"] = state_summary.fits / len(diagnostics)
    state_summary.to_csv(
        ROOT / "results" / "schafer_ten_pair_diagnostic_summary.csv", index=False
    )

    comparison_rows = []
    for name, column in {
        "current_tuned_l1_14_features": "supervised_auprc",
        "fixed_l2_mean_of_three_scores": "mean_three_scores_auprc",
        "fixed_l2_three_scores": "three_scores_auprc",
        "fixed_l2_three_scores_plus_joint_prevalence":
            "three_scores_plus_joint_prevalence_auprc",
        "fixed_unsupervised_mean_score": "fixed_mean_score_auprc",
        "adaptive_mean_score_calibration": "adaptive_simple_auprc",
        "OneNet": "onenet_auprc",
    }.items():
        for subset_name, subset in {
            "all fits": diagnostics,
            "at least one training positive": diagnostics[
                diagnostics.training_positives.gt(0)
            ],
        }.items():
            comparison_rows.append({
                "method": name,
                "subset": subset_name,
                "fits": len(subset),
                "mean_auprc": subset[column].mean(),
                "sd_auprc": subset[column].std(),
            })
    comparison = pd.DataFrame(comparison_rows)
    comparison.to_csv(
        ROOT / "results" / "schafer_ten_pair_complexity_comparison.csv", index=False
    )

    print(state_summary.to_string(index=False))
    print("\nSelected settings")
    print(diagnostics.selected_setting.value_counts().to_string())
    print("\nTraining-positive counts")
    print(diagnostics.training_positives.value_counts().sort_index().to_string())
    print("\nComplexity comparison")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
