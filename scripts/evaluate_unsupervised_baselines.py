#!/usr/bin/env python3
"""Evaluate unsupervised edge scores and native graph selections by test fold."""

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
METHODS = ("pln", "glmnet", "spieceasi", "spring", "sparcc")


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def average_precision(labels, scores):
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(zip(scores, labels), reverse=True)
    tp = fp = 0
    previous_recall = 0.0
    ap = 0.0
    index = 0
    while index < len(ordered):
        score = ordered[index][0]
        group_tp = group_fp = 0
        while index < len(ordered) and ordered[index][0] == score:
            if ordered[index][1] == 1:
                group_tp += 1
            else:
                group_fp += 1
            index += 1
        tp += group_tp
        fp += group_fp
        recall = tp / positives
        precision = tp / (tp + fp)
        ap += (recall - previous_recall) * precision
        previous_recall = recall
    return ap


def auroc(labels, scores):
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None
    ordered = sorted(zip(scores, labels))
    rank_sum = 0.0
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][0] == ordered[index][0]:
            end += 1
        average_rank = ((index + 1) + end) / 2
        rank_sum += average_rank * sum(label for _, label in ordered[index:end])
        index = end
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def native_metrics(labels, selected):
    tp = sum(y == 1 and s == 1 for y, s in zip(labels, selected))
    fp = sum(y == 0 and s == 1 for y, s in zip(labels, selected))
    fn = sum(y == 1 and s == 0 for y, s in zip(labels, selected))
    precision = tp / (tp + fp) if tp + fp else None
    recall = tp / (tp + fn) if tp + fn else None
    f1 = (2 * precision * recall / (precision + recall)
          if precision is not None and recall is not None and precision + recall else None)
    return tp + fp, precision, recall, f1


def display(value):
    return "NA" if value is None or not math.isfinite(value) else f"{value:.6f}"


def main():
    manifest = read_rows(ROOT / "training_data" / "fold_manifest.csv")
    results = []
    for fold in manifest:
        rows = read_rows(ROOT / "training_data" / fold["fold"] / "test.csv")
        labels = [int(row["interaction_label"]) for row in rows]
        for method in METHODS:
            scores = [float(row[f"{method}_score_percentile"]) for row in rows]
            selected = [int(float(row[f"{method}_selected"])) for row in rows]
            n_selected, precision, recall, f1 = native_metrics(labels, selected)
            results.append({
                "fold": fold["fold"],
                "held_out_studies": fold["held_out_studies"],
                "method": method,
                "n_test": len(labels),
                "n_positive": sum(labels),
                "n_neutral": len(labels) - sum(labels),
                "average_precision": display(average_precision(labels, scores)),
                "auroc": display(auroc(labels, scores)),
                "selected_edges": n_selected,
                "selected_density": display(n_selected / len(labels)),
                "precision": display(precision),
                "recall": display(recall),
                "f1": display(f1),
            })

    output_dir = ROOT / "results"
    output_dir.mkdir(exist_ok=True)
    output = output_dir / "unsupervised_baseline_by_fold.csv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)
    print(output)


if __name__ == "__main__":
    main()
