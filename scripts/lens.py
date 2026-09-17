"""LENS centred ridge model, abundance features, and training-only tuning."""
from __future__ import annotations
from pathlib import Path
from dataclasses import dataclass
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import average_precision_score

ROOT = Path(__file__).resolve().parents[1]
SEED = 20260821
SYSTEMS = ["butyrate_assembly_2021", "carlstrom_phyllosphere_2019", "schafer_phyllosphere_2022"]
# Preserve the seed offsets used in the saved manuscript experiments.
SYSTEM_SEED_INDEX = dict(zip(SYSTEMS, [0, 1, 3]))
BUDGET_FRACTIONS = [0.20, 0.40, 0.60, 0.80]
REFERENCE_METHODS = {s: ("OneNet", "onenet_score") for s in SYSTEMS}
FIXED_FALLBACK_FEATURES = ["pln_score_percentile", "glmnet_score_percentile", "sparcc_score_percentile"]


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

FEATURES = FIXED_FALLBACK_FEATURES + DIRECT_FEATURES

def canonicalize(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    first = frame["taxon_1"].astype(str)
    second = frame["taxon_2"].astype(str)
    frame["taxon_1"] = np.minimum(first, second)
    frame["taxon_2"] = np.maximum(first, second)
    return frame

def load_base(system):
    if system not in SYSTEMS:
        raise ValueError(f"Unsupported system: {system}")
    truth = canonicalize(pd.read_csv(ROOT / "cleaned_data" / f"{system}_tested_pairs.csv"))
    features = canonicalize(pd.read_csv(ROOT / "analysis_data" / f"{system}_pair_features.csv"))
    frame = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], how="inner", validate="many_to_one")
    frame = frame[frame.interaction_label.notna()].copy().reset_index(drop=True)
    frame["interaction_label"] = frame.interaction_label.astype(int)
    frame["pair_id"] = frame.taxon_1 + "||" + frame.taxon_2
    frame["mean_score_rank"] = frame[FIXED_FALLBACK_FEATURES].mean(axis=1)
    return frame


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
    frame = load_base(system)
    direct = build_direct_features(system)
    frame = frame.merge(direct, on=["taxon_1", "taxon_2"], how="left", validate="many_to_one")
    if frame[DIRECT_FEATURES].isna().all(axis=0).any():
        missing = frame[DIRECT_FEATURES].columns[frame[DIRECT_FEATURES].isna().all()].tolist()
        raise ValueError(f"Entire direct feature columns missing for {system}: {missing}")
    return frame

def load_system(system: str) -> pd.DataFrame:
    frame = load_analysis(system)
    onenet = canonicalize(pd.read_csv(
        ROOT / "analysis_data" / f"{system}_onenet_pair_scores.csv"
    ))
    frame = frame.merge(
        onenet[["dataset", "taxon_1", "taxon_2", "onenet_score"]],
        on=["dataset", "taxon_1", "taxon_2"],
        how="left",
        validate="many_to_one",
    )
    if frame["onenet_score"].isna().any():
        raise ValueError(f"{system}: missing OneNet scores after alignment")
    return frame

@dataclass(frozen=True)
class Candidate:
    family: str
    setting: str

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

    raise ValueError("LENS supports only centred ridge candidates")


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

CANDIDATES = [
    Candidate("centered_ridge", "lambda=0.01"),
    Candidate("centered_ridge", "lambda=0.1"),
    Candidate("centered_ridge", "lambda=1"),
    Candidate("centered_ridge", "lambda=10"),
]

def select_regularization(training: pd.DataFrame, seed: int):
    folds = grouped_folds(training, seed)
    if not folds:
        return CANDIDATES[0], []
    rows = []
    for candidate_index, candidate in enumerate(CANDIDATES):
        scores = []
        for fold_index, (inner_train, inner_valid) in enumerate(folds):
            y = inner_valid.interaction_label.to_numpy(int)
            if np.unique(y).size < 2:
                continue
            predicted = fit_score(
                candidate, inner_train, inner_valid,
                seed + 1000 * candidate_index + fold_index,
            )
            scores.append(average_precision_score(y, predicted))
        rows.append({
            "setting": candidate.setting,
            "inner_folds": len(scores),
            "mean_inner_auprc": float(np.mean(scores)) if scores else np.nan,
        })
    valid = [row for row in rows if np.isfinite(row["mean_inner_auprc"])]
    if not valid:
        return CANDIDATES[0], rows
    # Prefer stronger regularization when inner AUPRC ties.
    best = max(valid, key=lambda row: (row["mean_inner_auprc"],
                                       float(row["setting"].split("=")[1])))
    return Candidate("centered_ridge", best["setting"]), rows
