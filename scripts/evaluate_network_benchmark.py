#!/usr/bin/env python3
"""Evaluate unsupervised networks on one cleaned experimental benchmark."""

from pathlib import Path
import csv
import sys

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
DATASET = sys.argv[1] if len(sys.argv) > 1 else "wortel_syncom_2026"
METHODS = ["pln", "glmnet", "spieceasi", "spring", "sparcc"]


def classification_metrics(labels, selected):
    labels = pd.Series(labels).astype(int)
    selected = pd.Series(selected).astype(int)
    tp = int(((labels == 1) & (selected == 1)).sum())
    fp = int(((labels == 0) & (selected == 1)).sum())
    fn = int(((labels == 1) & (selected == 0)).sum())
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else float("nan")
    return int(selected.sum()), precision, recall, f1


def main():
    features = pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_pair_features.csv")
    truth = pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_tested_pairs.csv")
    for frame in (features, truth):
        ordered = frame[["taxon_1", "taxon_2"]].apply(lambda row: sorted(row.astype(str)), axis=1, result_type="expand")
        frame[["taxon_1", "taxon_2"]] = ordered.to_numpy()
    status = pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_method_status.csv").set_index("method")
    data = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one")
    labels = data["interaction_label"].astype(int)
    edge_budget = int(labels.sum())
    rows = []

    abundance = pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_abundance.csv")
    for method in METHODS:
        scores = data[f"{method}_score_percentile"].astype(float)
        native = data[f"{method}_selected"].astype(int)
        native_n, precision, recall, f1 = classification_metrics(labels, native)

        ranked = data.assign(_score=scores).sort_values(
            ["_score", "taxon_1", "taxon_2"], ascending=[False, True, True]
        )
        density_selected = pd.Series(0, index=data.index)
        density_selected.loc[ranked.index[:edge_budget]] = 1
        _, density_precision, density_recall, density_f1 = classification_metrics(labels, density_selected)

        rows.append({
            "method": method,
            "n_samples": len(abundance),
            "n_taxa": abundance.shape[1] - 1,
            "tested_pairs": len(data),
            "positive_pairs": edge_budget,
            "neutral_pairs": int((labels == 0).sum()),
            "average_precision": average_precision_score(labels, scores),
            "auroc": roc_auc_score(labels, scores),
            "native_selected_edges": native_n,
            "native_density": native_n / len(data),
            "native_precision": precision,
            "native_recall": recall,
            "native_f1": f1,
            "density_matched_edges": edge_budget,
            "density_matched_precision": density_precision,
            "density_matched_recall": density_recall,
            "density_matched_f1": density_f1,
            "runtime_seconds": status.loc[method, "elapsed_seconds"],
        })

    output = ROOT / "results" / f"{DATASET}_network_performance.csv"
    pd.DataFrame(rows).to_csv(output, index=False, quoting=csv.QUOTE_MINIMAL)
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
