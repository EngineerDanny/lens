#!/usr/bin/env python3
"""Independently verify saved split rules and discrimination metrics."""
import json
import numpy as np
import pandas as pd
from lens import ROOT, SYSTEMS

OUT = ROOT / "results/measurement_placement"


def independent_ap(y, s):
    order = np.argsort(-s, kind="stable")
    y, s = y[order], s[order]
    ends = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s) - 1]
    tp = np.cumsum(y)[ends]
    recall = tp / y.sum()
    precision = tp / (ends + 1)
    return np.sum(np.diff(np.r_[0, recall]) * precision)


def independent_auc(y, s):
    ranks = pd.Series(s).rank(method="average").to_numpy()
    p, n = y.sum(), (1 - y).sum()
    return (ranks[y == 1].sum() - p * (p + 1) / 2) / (p * n)


def main():
    metrics = pd.read_csv(OUT / "partition_metrics.csv")
    cols = {
        "LENS restricted": "lens_restricted", "LENS distributed": "lens_distributed",
        "OneNet": "onenet_score", "Fixed mean": "reference_mean",
        "PLNNetwork": "pln_score_percentile", "Poisson GLMNet": "glmnet_score_percentile",
        "SparCC": "sparcc_score_percentile", "SPIEC-EASI": "spieceasi_score_percentile",
    }
    checked = 0
    for system in SYSTEMS:
        plans = json.loads((OUT / f"{system}_plan.json").read_text())
        by_repeat = {p["repeat"]: p for p in plans["partitions"]}
        predictions = pd.read_csv(OUT / f"{system}_predictions.csv")
        subset = metrics[metrics.analysis_set == system]
        for (repeat, budget), group in predictions.groupby(["repeat", "budget_pairs"]):
            plan = by_repeat[repeat]
            assert set(group.pair_id) == set(plan["test_ids"])
            assert group.pair_id.is_unique
            reserve = set(plan["reserved_taxa"])
            for design in ("restricted", "distributed"):
                training = plan[f"{design}_order"][:budget]
                assert len(set(training)) == budget
                assert not set(training) & set(group.pair_id)
                if design == "restricted":
                    assert not {t for pair in training for t in pair.split("||")} & reserve
            y = group.interaction_label.to_numpy(int)
            for method, column in cols.items():
                row = subset[(subset.repeat == repeat) & (subset.budget_pairs == budget) & (subset.method == method)].iloc[0]
                scores = group[column].to_numpy(float)
                np.testing.assert_allclose(independent_ap(y, scores), row.auprc, atol=1e-12)
                np.testing.assert_allclose(independent_auc(y, scores), row.auroc, atol=1e-12)
                checked += 1
            if budget == 5:
                np.testing.assert_allclose(group.lens_restricted, group.reference_mean)
                np.testing.assert_allclose(group.lens_distributed, group.reference_mean)
        static = subset[subset.method.isin(["OneNet", "Fixed mean"])]
        assert (static.groupby(["repeat", "method"]).auprc.nunique() == 1).all()
        for method in ["LENS distributed", "LENS restricted", "OneNet"]:
            sets = subset[subset.method == method].groupby("budget_pairs")["repeat"].apply(set).tolist()
            assert all(s == sets[0] for s in sets)
        print(f"PASS {system}: split isolation, identical tests across budgets, fixed references", flush=True)
    print(f"PASS {checked} method evaluations independently recomputed for AP and AUROC")


if __name__ == "__main__":
    main()
