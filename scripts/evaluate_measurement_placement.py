#!/usr/bin/env python3
"""Compare restricted and distributed labels on identical hidden taxon pairs.

The experiment leaves scripts/lens.py and all existing results unchanged.
Partitions and budget caps use taxon identifiers and pair availability only.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

from lens import ROOT, SEED, SYSTEMS, SYSTEM_SEED_INDEX
from lens import FIXED_FALLBACK_FEATURES, load_system, select_regularization, fit_score

OUT = ROOT / "results" / "measurement_placement"
REPEATS = 20
RESERVED_FRACTION = 0.4
GRID = [5, 10, 20, 50, 100, 200, 500, 1000]
BASELINES = {
    "OneNet": "onenet_score",
    "Fixed mean": "reference_mean",
    "PLNNetwork": "pln_score_percentile",
    "Poisson GLMNet": "glmnet_score_percentile",
    "SparCC": "sparcc_score_percentile",
    "SPIEC-EASI": "spieceasi_score_percentile",
}


def plan_partitions(frame, system):
    """Generate plans without consulting any experimental outcomes."""
    taxa = np.array(sorted(set(frame.taxon_1) | set(frame.taxon_2)))
    plans = []
    for repeat in range(REPEATS):
        seed = SEED + 70_000_000 + 10_000 * SYSTEM_SEED_INDEX[system] + repeat
        rng = np.random.default_rng(seed)
        reserved = set(rng.choice(taxa, max(2, round(RESERVED_FRACTION * len(taxa))), replace=False))
        first, second = frame.taxon_1.isin(reserved), frame.taxon_2.isin(reserved)
        test = frame.loc[first & second, "pair_id"].tolist()
        restricted = set(frame.loc[~first & ~second, "pair_id"])
        # Shared random priorities couple the designs; each eligible ordering
        # remains a uniform random permutation of its own pool.
        order = rng.permutation(sorted(set(frame.pair_id) - set(test))).tolist()
        restricted_order = [pair for pair in order if pair in restricted]
        plans.append({
            "repeat": repeat + 1, "seed": seed, "reserved_taxa": sorted(reserved),
            "test_ids": test, "distributed_order": order,
            "restricted_order": restricted_order,
        })
    cap = min(len(p["restricted_order"]) for p in plans)
    if cap < 1:
        raise ValueError(f"{system}: no positive common label budget across partitions")
    # The cap is shared by every partition, so increasing n never changes
    # the collection of test sets included in a system's curve.
    budgets = sorted(set([n for n in GRID if n <= cap] + [cap]))
    return plans, budgets


def metrics(labels, scores):
    return {
        "auprc": average_precision_score(labels, scores),
        "auroc": roc_auc_score(labels, scores),
    }


def run_system(system):
    frame = load_system(system)
    assert frame.pair_id.is_unique
    assert frame.interaction_label.isin([0, 1]).all()
    frame["reference_mean"] = frame[FIXED_FALLBACK_FEATURES].mean(axis=1)
    plans, budgets = plan_partitions(frame, system)
    # Guard against accidental use of label values in partition construction.
    altered = frame.copy()
    altered["interaction_label"] = 1 - altered.interaction_label
    assert plan_partitions(altered, system) == (plans, budgets)
    (OUT / f"{system}_plan.json").write_text(json.dumps({
        "system": system, "budgets": budgets, "partitions": plans,
    }, indent=2))
    records, audit, predictions, tuning = [], [], [], []
    for plan in plans:
        repeat = plan["repeat"]
        test = frame[frame.pair_id.isin(plan["test_ids"])].copy()
        y = test.interaction_label.to_numpy(int)
        status = "ok" if len(np.unique(y)) == 2 else (
            "no_test_pairs" if not len(y) else "one_test_class")
        common = {
            "analysis_set": system, "repeat": repeat, "status": status,
            "test_pairs": len(test), "test_interactions": int(y.sum()),
            "test_neutrals": int((y == 0).sum()),
        }
        audit.append({**common, "reserved_taxa": len(plan["reserved_taxa"]),
                      "test_taxa_with_pairs": len(set(test.taxon_1) | set(test.taxon_2)),
                      "restricted_pool": len(plan["restricted_order"]),
                      "distributed_pool": len(plan["distributed_order"])})
        if status != "ok":
            print(system, repeat, status, flush=True)
            continue
        static = {name: metrics(y, test[col].to_numpy(float)) for name, col in BASELINES.items()}
        for budget in budgets:
            for name, value in static.items():
                records.append({**common, "budget_pairs": budget, "method": name, **value})
            pair_predictions = test[["pair_id", "taxon_1", "taxon_2", "interaction_label", *BASELINES.values()]].copy()
            pair_predictions["analysis_set"] = system
            pair_predictions["repeat"] = repeat
            pair_predictions["budget_pairs"] = budget
            for design in ["restricted", "distributed"]:
                ids = plan[f"{design}_order"][:budget]
                training = frame[frame.pair_id.isin(ids)].copy()
                assert len(training) == budget
                assert not set(training.pair_id) & set(test.pair_id)
                covered = set(training.taxon_1) | set(training.taxon_2)
                if design == "restricted":
                    assert not covered & set(plan["reserved_taxa"])
                counts = training.interaction_label.value_counts().reindex([0, 1], fill_value=0)
                fallback = counts.min() < 3
                tuning_seed = plan["seed"] + 100_000 + budget
                try:
                    candidate, inner = select_regularization(training, tuning_seed)
                    score = fit_score(candidate, training, test, tuning_seed)
                except Exception as error:
                    # A fitting failure must not disappear from the summary.
                    raise RuntimeError(f"{system} repeat {repeat} n={budget} {design}: {error}") from error
                if fallback:
                    np.testing.assert_allclose(score, test.reference_mean)
                for row in inner:
                    tuning.append({"analysis_set": system, "repeat": repeat,
                                   "budget_pairs": budget, "design": design, **row})
                test_taxa = set(test.taxon_1) | set(test.taxon_2)
                records.append({
                    **common, "budget_pairs": budget, "method": f"LENS {design}",
                    **metrics(y, score), "fallback": bool(fallback),
                    "selected_setting": candidate.setting,
                    "training_interactions": int(counts.loc[1]),
                    "training_neutrals": int(counts.loc[0]),
                    "test_taxon_coverage": len(test_taxa & covered) / len(test_taxa),
                    "test_pairs_both_taxa_covered": float((test.taxon_1.isin(covered) & test.taxon_2.isin(covered)).mean()),
                })
                pair_predictions[f"lens_{design}"] = score
            predictions.extend(pair_predictions.to_dict("records"))
        print(system, f"{repeat}/{REPEATS}", f"test={len(test)}, interactions={y.sum()}, budgets={budgets}", flush=True)
        # Save progress after each complete partition; never reuse stale fits.
        pd.DataFrame(records).to_csv(OUT / f"{system}_metrics.csv", index=False)
    pd.DataFrame(audit).to_csv(OUT / f"{system}_audit.csv", index=False)
    pd.DataFrame(predictions).to_csv(OUT / f"{system}_predictions.csv", index=False)
    pd.DataFrame(tuning).to_csv(OUT / f"{system}_tuning.csv", index=False)
    return records, audit


def summarize(records):
    frame = pd.DataFrame(records)
    summary = frame.groupby(["analysis_set", "budget_pairs", "method"]).agg(
        partitions=("repeat", "size"), mean_auprc=("auprc", "mean"),
        sd_auprc=("auprc", "std"), mean_auroc=("auroc", "mean"),
        sd_auroc=("auroc", "std"), mean_test_pairs=("test_pairs", "mean"),
        fallback_rate=("fallback", "mean"),
        mean_training_interactions=("training_interactions", "mean"),
        mean_test_taxon_coverage=("test_taxon_coverage", "mean"),
    ).reset_index()
    frame.to_csv(OUT / "partition_metrics.csv", index=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    largest = summary[summary.budget_pairs == summary.groupby("analysis_set").budget_pairs.transform("max")]
    largest.to_csv(OUT / "largest_budget_summary.csv", index=False)
    comparisons = [("LENS distributed", "LENS restricted"),
                   ("LENS restricted", "Fixed mean"),
                   ("LENS distributed", "Fixed mean"),
                   ("LENS restricted", "OneNet"), ("LENS distributed", "OneNet")]
    differences = []
    for metric in ["auprc", "auroc"]:
        pivot = frame.pivot(index=["analysis_set", "budget_pairs", "repeat"], columns="method", values=metric)
        for first, second in comparisons:
            delta = (pivot[first] - pivot[second]).rename("difference").reset_index()
            delta["comparison"] = f"{first} minus {second}"
            delta["metric"] = metric
            differences.append(delta)
    paired = pd.concat(differences, ignore_index=True)
    paired.to_csv(OUT / "paired_differences.csv", index=False)
    paired.groupby(["analysis_set", "budget_pairs", "comparison", "metric"]).agg(
        partitions=("repeat", "size"), mean_difference=("difference", "mean"),
        sd_difference=("difference", "std"), median_difference=("difference", "median"),
        win_fraction=("difference", lambda x: np.mean(x > 1e-12)),
        tie_fraction=("difference", lambda x: np.mean(np.abs(x) <= 1e-12)),
    ).reset_index().to_csv(OUT / "paired_summary.csv", index=False)
    return summary


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [ROOT / "scripts/lens.py"]
    for system in SYSTEMS:
        inputs += [ROOT / "cleaned_data" / f"{system}_{kind}.csv" for kind in ["abundance", "tested_pairs"]]
        inputs += [ROOT / "analysis_data" / f"{system}_{kind}.csv" for kind in ["pair_features", "onenet_pair_scores"]]
    metadata = {
        "seed": SEED, "repeats": REPEATS, "reserved_taxon_fraction": RESERVED_FRACTION,
        "test": "all observed tested pairs between two reserved taxa",
        "restricted": "training pairs between two nonreserved taxa",
        "distributed": "all remaining tested pairs; exact test pairs excluded",
        "sampling": "shared uniform random priorities; nested; labels never used",
        "budgets": "grid capped at minimum restricted pool size across all predetermined partitions",
        "invalid_tests": "recorded without redrawing; no metric if fewer than two test classes",
        "model": "unchanged centred ridge LENS with existing inner pair CV and fallback",
        "interval": "boxplots across eligible taxon partitions: median, quartiles, whiskers within 1.5 IQR, and outliers; not confidence intervals",
        "scope": "within-system taxa without interaction labels, but with abundance observations",
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2))
    records, audit = [], []
    for system in SYSTEMS:
        rows, checks = run_system(system)
        records.extend(rows)
        audit.extend(checks)
    pd.DataFrame(audit).to_csv(OUT / "partition_audit.csv", index=False)
    summary = summarize(records)
    print(summary[summary.method.isin(["LENS restricted", "LENS distributed", "OneNet", "Fixed mean"])].to_string(index=False))


if __name__ == "__main__":
    main()
