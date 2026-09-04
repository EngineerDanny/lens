#!/usr/bin/env python3
"""Parse and evaluate the CoNet result on its prespecified taxon universe."""

from pathlib import Path
import csv
import math

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
DATASET = "schafer_phyllosphere_2022"
OTHER_METHODS = ["pln", "glmnet", "spieceasi", "spring", "sparcc"]


def canonical_pair(a, b):
    return tuple(sorted((str(a), str(b))))


def metrics(labels, selected):
    labels = pd.Series(labels).astype(int)
    selected = pd.Series(selected).astype(int)
    tp = int(((labels == 1) & (selected == 1)).sum())
    fp = int(((labels == 0) & (selected == 1)).sum())
    fn = int(((labels == 1) & (selected == 0)).sum())
    precision = tp / (tp + fp) if tp + fp else math.nan
    recall = tp / (tp + fn) if tp + fn else math.nan
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else math.nan
    return int(selected.sum()), precision, recall, f1


def add_result(rows, method, data, scores, selected, runtime):
    labels = data["interaction_label"].astype(int)
    edge_budget = int(labels.sum())
    native_n, precision, recall, f1 = metrics(labels, selected)
    order = pd.DataFrame({"score": scores, "a": data.taxon_1, "b": data.taxon_2}).sort_values(
        ["score", "a", "b"], ascending=[False, True, True]
    )
    density = pd.Series(0, index=data.index)
    density.loc[order.index[:edge_budget]] = 1
    _, dp, dr, df = metrics(labels, density)
    rows.append({
        "method": method,
        "n_samples": 662,
        "n_taxa": 53,
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
        "density_matched_precision": dp,
        "density_matched_recall": dr,
        "density_matched_f1": df,
        "runtime_seconds": runtime,
    })


def main():
    abundance = pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_abundance.csv")
    eligible = set(abundance.columns[1:][(abundance.iloc[:, 1:] > 0).sum(axis=0) >= 20])

    truth = pd.read_csv(ROOT / "cleaned_data" / f"{DATASET}_tested_pairs.csv")
    truth[["taxon_1", "taxon_2"]] = truth.apply(
        lambda row: canonical_pair(row.taxon_1, row.taxon_2), axis=1, result_type="expand"
    )
    truth = truth[truth.taxon_1.isin(eligible) & truth.taxon_2.isin(eligible)].copy()

    features = pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_pair_features.csv")
    features[["taxon_1", "taxon_2"]] = features.apply(
        lambda row: canonical_pair(row.taxon_1, row.taxon_2), axis=1, result_type="expand"
    )
    data = truth.merge(features, on=["dataset", "taxon_1", "taxon_2"], validate="one_to_one")

    conet = pd.read_csv(
        ROOT / "analysis_cache" / DATASET / "conet" / "conet_minocc20_100.tab",
        sep="\t", comment="#"
    )
    endpoints = conet.Label.str.split("->", n=1, expand=True)
    conet[["taxon_1", "taxon_2"]] = endpoints.apply(
        lambda row: canonical_pair(row.iloc[0], row.iloc[1]), axis=1, result_type="expand"
    )
    conet = conet.drop_duplicates(["taxon_1", "taxon_2"])
    conet["conet_selected"] = 1
    # CoNet's 100-resample p-values underflow for most retained edges. Method
    # support supplies the primary ordering; the corrected p-value breaks ties.
    q = conet.qval.clip(lower=1e-300)
    conet["conet_score"] = conet.method_number.astype(float) + (-q.map(math.log10) / 1000)
    data = data.merge(
        conet[["taxon_1", "taxon_2", "conet_selected", "conet_score", "qval", "method_number"]],
        on=["taxon_1", "taxon_2"], how="left", validate="one_to_one"
    )
    data[["conet_selected", "conet_score", "method_number"]] = data[
        ["conet_selected", "conet_score", "method_number"]
    ].fillna(0)
    data["conet_qvalue"] = data.qval

    status = pd.read_csv(ROOT / "analysis_data" / f"{DATASET}_method_status.csv").set_index("method")
    rows = []
    for method in OTHER_METHODS:
        add_result(
            rows, method, data,
            data[f"{method}_score_percentile"].astype(float),
            data[f"{method}_selected"].astype(int),
            float(status.loc[method, "elapsed_seconds"]),
        )
    add_result(rows, "CoNet", data, data.conet_score, data.conet_selected, 138.2)

    columns = [
        "dataset", "taxon_1", "taxon_2", "interaction_label", "tested_status",
        "conet_selected", "conet_score", "conet_qvalue", "method_number"
    ]
    data[columns].to_csv(ROOT / "analysis_data" / f"{DATASET}_conet_tested_pair_scores.csv", index=False)
    result = pd.DataFrame(rows)
    result.to_csv(ROOT / "results" / f"{DATASET}_network_performance_conet_universe.csv", index=False,
                  quoting=csv.QUOTE_MINIMAL)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
