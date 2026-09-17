#!/usr/bin/env python3
"""Exploratory density recovery on the existing Carlstrom 80% outer folds.

LENS thresholds maximize F1 on nested training predictions, with higher
thresholds breaking ties. OneNet selects mean frequencies strictly above 0.9
without interaction labels.
Density is defined only among experimentally tested pairs.
"""
import json
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

from evaluate_pair_five_fold_cv import load_system
from evaluate_optimized_sparse_logistic import select_regularization
from evaluate_within_system_calibration import ROOT, SEED
from optimize_supervised_model import fit_score, grouped_folds


def threshold_from_training(y, score):
    assert np.isfinite(score).all()
    precision, recall, thresholds = precision_recall_curve(y, score)
    f1 = 2 * precision[:-1] * recall[:-1] / np.maximum(
        precision[:-1] + recall[:-1], 1e-15)
    best = np.flatnonzero(np.isclose(f1, f1.max(), rtol=0, atol=1e-12))
    return float(thresholds[best[-1]])


def main():
    system = "carlstrom_phyllosphere_2019"
    frame = load_system(system)
    saved = pd.read_csv(ROOT / "results/final_sparse_pr_predictions_80pct.csv")
    saved = saved[saved.analysis_set == system].copy()
    assert frame.pair_id.is_unique and saved.pair_id.is_unique
    assert set(frame.pair_id) == set(saved.pair_id)
    folds, predictions, inner_records = [], [], []
    for fold in sorted(saved.fold.unique()):
        test = saved[saved.fold == fold].copy()
        training = frame[~frame.pair_id.isin(test.pair_id)].copy()
        assert set(training.pair_id).isdisjoint(test.pair_id)
        assert (test.budget_pairs == len(training)).all()
        assert np.array_equal(test.interaction_label.to_numpy(),
                              frame.set_index("pair_id").loc[test.pair_id, "interaction_label"].to_numpy())
        inner = []
        for index, (fit, valid) in enumerate(grouped_folds(training, SEED + 900_000 + fold)):
            seed = SEED + 910_000 + fold * 100 + index
            candidate, _ = select_regularization(fit, seed)
            score = fit_score(candidate, fit, valid, seed)
            record = valid[["pair_id", "interaction_label"]].copy()
            record["supervised_score"] = score
            record["outer_fold"] = fold
            record["inner_fold"] = index + 1
            record["selected_setting"] = candidate.setting
            inner.append(record)
        calibration = pd.concat(inner, ignore_index=True)
        assert calibration.pair_id.is_unique
        assert set(calibration.pair_id) == set(training.pair_id)
        inner_records.append(calibration)
        y = test.interaction_label.to_numpy(int)
        baseline = float(training.interaction_label.mean())
        for method, column in [("LENS", "supervised_score"), ("OneNet", "onenet_score")]:
            if method == "LENS":
                threshold = threshold_from_training(calibration.interaction_label, calibration[column])
                selected = test[column].to_numpy() >= threshold
                comparison = ">="
            else:
                threshold = 0.9
                selected = test[column].to_numpy() > threshold
                comparison = ">"
            tp = int(np.sum(selected & (y == 1)))
            fp = int(np.sum(selected & (y == 0)))
            fn = int(np.sum(~selected & (y == 1)))
            folds.append(dict(fold=int(fold), method=method, test_pairs=len(y),
                true_interactions=int(y.sum()), selected_interactions=int(selected.sum()),
                true_density=float(y.mean()), predicted_density=float(selected.mean()),
                absolute_density_error=float(abs(selected.mean() - y.mean())),
                threshold=threshold, comparison=comparison, tp=tp, fp=fp, fn=fn,
                training_frequency=baseline,
                training_frequency_absolute_error=float(abs(baseline - y.mean()))))
            out = test[["pair_id", "fold", "interaction_label"]].copy()
            out["method"] = method
            out["selected"] = selected
            out["score"] = test[column].to_numpy()
            out["threshold"] = threshold
            out["comparison"] = comparison
            predictions.append(out)
        print(f"Finished outer fold {fold}", flush=True)
    result = pd.DataFrame(folds)
    summary = []
    for method, group in result.groupby("method", sort=False):
        tp, fp, fn = group[["tp", "fp", "fn"]].sum()
        summary.append(dict(method=method, tested_pairs=int(group.test_pairs.sum()),
            true_interactions=int(group.true_interactions.sum()),
            selected_interactions=int(group.selected_interactions.sum()),
            predicted_density=group.selected_interactions.sum()/group.test_pairs.sum(),
            absolute_pooled_density_error=abs(group.selected_interactions.sum()-group.true_interactions.sum())/group.test_pairs.sum(),
            mean_fold_absolute_density_error=group.absolute_density_error.mean(),
            precision=tp/(tp+fp) if tp+fp else np.nan,
            recall=tp/(tp+fn), f1=2*tp/(2*tp+fp+fn)))
    output = ROOT / "results/carlstrom_density_recovery"
    output.mkdir(exist_ok=True)
    result.to_csv(output / "folds.csv", index=False)
    pd.concat(predictions).to_csv(output / "predictions.csv", index=False)
    pd.concat(inner_records).to_csv(output / "threshold_training_predictions.csv", index=False)
    pd.DataFrame(summary).to_csv(output / "summary.csv", index=False)
    metadata = dict(system=system, budget=0.8, seed=SEED,
        lens_threshold_rule="Maximum nested training F1; ties use highest threshold; select scores >= threshold",
        onenet_threshold_rule="Mean selection frequency > 0.9; no interaction labels",
        onenet_alignment_target=0.8,
        onenet_resamples=30,
        onenet_estimators=6,
        selection_reference="https://doi.org/10.1371/journal.pcbi.1012627",
        baseline_mean_absolute_density_error=float(result.drop_duplicates("fold").training_frequency_absolute_error.mean()),
        scope="989 experimentally tested pairs; no inference of truth for untested pairs",
        qualification="Existing label-stratified outer folds constrain prevalence differences. This is an exploratory diagnostic, not independent validation of density under distribution shift.",
        onenet="Network scores and selection threshold use no experimental interaction labels.",
        undefined_precision="Missing when no pairs are selected; recall and F1 remain zero.")
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
