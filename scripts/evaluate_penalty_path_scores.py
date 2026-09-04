#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score


ROOT = Path(__file__).resolve().parents[1]
DATASETS = [
    "butyrate_assembly_2021",
    "carlstrom_phyllosphere_2019",
    "schafer_phyllosphere_2022",
    "friedman_microcosm_2017",
]
METHODS = ["pln", "spieceasi", "spring"]
PATH_SUFFIXES = [
    "entry",
    "presence",
    "mean_abs",
    "max_abs",
    "terminal_abs",
    "stability_mean",
    "stability_max",
]


def canonical(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    ordered = np.sort(frame[["taxon_1", "taxon_2"]].astype(str).to_numpy(), axis=1)
    frame["taxon_1"] = ordered[:, 0]
    frame["taxon_2"] = ordered[:, 1]
    return frame


audit_rows = []
performance_rows = []

for dataset in DATASETS:
    features = canonical(pd.read_csv(ROOT / "analysis_data" / f"{dataset}_pair_features.csv"))
    path_scores = canonical(pd.read_csv(ROOT / "analysis_data" / f"{dataset}_penalty_path_scores.csv"))
    truth_name = (
        f"{dataset}_tested_pairs_directional.csv"
        if dataset == "friedman_microcosm_2017"
        else f"{dataset}_tested_pairs.csv"
    )
    truth = canonical(pd.read_csv(ROOT / "cleaned_data" / truth_name))
    if dataset != "friedman_microcosm_2017":
        truth = truth.drop_duplicates(["taxon_1", "taxon_2"])
    data = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], how="inner")
    data = data.merge(path_scores, on=["dataset", "taxon_1", "taxon_2"], how="inner")
    y = data["interaction_label"].astype(int).to_numpy()

    for method in ["pln", "glmnet", "spieceasi", "spring", "sparcc"]:
        score_col = f"{method}_abs_score"
        score = data[score_col].fillna(0).astype(float).to_numpy()
        audit_rows.append({
            "dataset": dataset,
            "method": method,
            "tested_pairs": len(data),
            "positive_pairs": int(y.sum()),
            "zero_score_pairs": int(np.isclose(score, 0).sum()),
            "zero_score_fraction": float(np.isclose(score, 0).mean()),
            "unique_score_values": int(pd.Series(score).nunique()),
            "largest_tie_fraction": float(pd.Series(score).value_counts(normalize=True).iloc[0]),
            "current_auprc": float(average_precision_score(y, score)) if len(np.unique(y)) > 1 else np.nan,
        })

    if len(np.unique(y)) < 2:
        continue
    for method in METHODS:
        current = data[f"{method}_abs_score"].fillna(0).astype(float).to_numpy()
        performance_rows.append({
            "dataset": dataset,
            "method": method,
            "score": "selected_graph_weight",
            "tested_pairs": len(data),
            "positive_pairs": int(y.sum()),
            "random_auprc": float(y.mean()),
            "average_precision": float(average_precision_score(y, current)),
        })
        for suffix in PATH_SUFFIXES:
            col = f"{method}_path_{suffix}"
            if col not in data:
                continue
            score = data[col].fillna(0).astype(float).to_numpy()
            performance_rows.append({
                "dataset": dataset,
                "method": method,
                "score": f"path_{suffix}",
                "tested_pairs": len(data),
                "positive_pairs": int(y.sum()),
                "random_auprc": float(y.mean()),
                "average_precision": float(average_precision_score(y, score)),
            })

    def percentile(values: pd.Series) -> np.ndarray:
        return values.astype(float).rank(method="average", pct=True).to_numpy()

    current_consensus = np.mean([
        percentile(data[f"{method}_abs_score"].fillna(0))
        for method in ["pln", "spieceasi", "spring", "sparcc"]
    ], axis=0)
    path_consensus = np.mean([
        percentile(data["pln_path_mean_abs"].fillna(0)),
        percentile(data["spieceasi_path_mean_abs"].fillna(0)),
        percentile(data["spring_path_mean_abs"].fillna(0)),
        percentile(data["sparcc_abs_score"].fillna(0)),
    ], axis=0)
    for name, score in [
        ("current_four_method_mean_rank", current_consensus),
        ("path_integrated_four_method_mean_rank", path_consensus),
    ]:
        performance_rows.append({
            "dataset": dataset,
            "method": "consensus",
            "score": name,
            "tested_pairs": len(data),
            "positive_pairs": int(y.sum()),
            "random_auprc": float(y.mean()),
            "average_precision": float(average_precision_score(y, score)),
        })

audit = pd.DataFrame(audit_rows)
performance = pd.DataFrame(performance_rows)
audit.to_csv(ROOT / "results" / "continuous_score_sparsity_audit.csv", index=False)
performance.to_csv(ROOT / "results" / "penalty_path_score_performance.csv", index=False)

best = (
    performance.sort_values("average_precision", ascending=False)
    .groupby(["dataset", "method"], as_index=False)
    .first()
)
best.to_csv(ROOT / "results" / "penalty_path_score_best_diagnostic.csv", index=False)

summary_scores = [
    "selected_graph_weight",
    "path_mean_abs",
    "current_four_method_mean_rank",
    "path_integrated_four_method_mean_rank",
]
summary = performance[performance["score"].isin(summary_scores)].copy()
summary.to_csv(ROOT / "results" / "penalty_path_score_summary.csv", index=False)

print("\nSparsity audit")
print(audit.to_string(index=False))
print("\nBest diagnostic score per method and system")
print(best.to_string(index=False))
