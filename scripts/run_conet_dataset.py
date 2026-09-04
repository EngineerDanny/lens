#!/usr/bin/env python3
"""Run CoNet for one cleaned benchmark dataset and evaluate tested pairs."""

from __future__ import annotations

import argparse
import math
import subprocess
import time
from pathlib import Path

import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score


ROOT = Path(__file__).resolve().parents[1]
JAR = ROOT / "external_data" / "conet" / "conet-1.1.1.beta.jar"
MAIN_CLASS = "be.ac.vub.bsb.cooccurrence.cmd.CooccurrenceAnalyser"
METHODS = "correl_pearson/correl_spearman/dist_bray/dist_kullbackleibler/sim_mutInfo"


def canonical_pair(a: object, b: object) -> tuple[str, str]:
    return tuple(sorted((str(a), str(b))))


def run_command(args: list[str]) -> None:
    subprocess.run(args, check=True)


def prepare_input(dataset: str, cache_dir: Path, min_occurrence: int) -> tuple[int, int, int]:
    abundance_path = ROOT / "cleaned_data" / f"{dataset}_abundance.csv"
    abundance = pd.read_csv(abundance_path)
    occurrences = (abundance.iloc[:, 1:] > 0).sum(axis=0)
    eligible_taxa = int((occurrences >= min_occurrence).sum())
    possible_pairs = eligible_taxa * (eligible_taxa - 1) // 2
    output_path = cache_dir / f"abundance_minocc{min_occurrence}.tsv"
    run_command(
        [
            "python3",
            str(ROOT / "scripts" / "prepare_conet_input.py"),
            str(abundance_path),
            str(output_path),
            "--min-occurrence",
            str(min_occurrence),
        ]
    )
    return int(abundance.shape[0]), eligible_taxa, possible_pairs


def run_conet(
    dataset: str,
    cache_dir: Path,
    min_occurrence: int,
    iterations: int,
    max_edges: int,
    force: bool,
) -> float:
    input_path = cache_dir / f"abundance_minocc{min_occurrence}.tsv"
    threshold_path = cache_dir / f"thresholds_minocc{min_occurrence}.txt"
    permutation_path = cache_dir / f"permutation_minocc{min_occurrence}_{iterations}.txt"
    permutation_gml = cache_dir / f"permutation_minocc{min_occurrence}_{iterations}.gml"
    bootstrap_path = cache_dir / f"bootstrap_minocc{min_occurrence}_{iterations}.txt"
    output_prefix = cache_dir / f"conet_minocc{min_occurrence}_{iterations}"
    tab_path = output_prefix.with_suffix(".tab")

    if tab_path.exists() and not force:
        return math.nan

    start = time.perf_counter()
    java_prefix = ["java", "-Xmx8g", "-cp", str(JAR), MAIN_CLASS]
    run_command(
        java_prefix
        + [
            "--input",
            str(input_path),
            "--matrixtype",
            "count",
            "--method",
            "ensemble",
            "--ensemblemethods",
            METHODS,
            "--thresholdguessing",
            "edgeNumber",
            "--guessingparam",
            str(max_edges),
            "--topbottom",
            "--output",
            str(threshold_path),
            "--verbosity",
            "fatal",
        ]
    )
    run_command(
        java_prefix
        + [
            "--input",
            str(input_path),
            "--matrixtype",
            "count",
            "--method",
            "ensemble",
            "--ensemblemethods",
            METHODS,
            "--ensembleparamfile",
            str(threshold_path),
            "--networkmergestrategy",
            "union",
            "--multigraph",
            "--filter",
            "rand",
            "--randroutine",
            "edgeScores",
            "--resamplemethod",
            "shuffle_rows",
            "--renorm",
            "--iterations",
            str(iterations),
            "--edgethreshold",
            "0.05",
            "--randscorefile",
            str(permutation_path),
            "--scoreexport",
            "--format",
            "gml",
            "--output",
            str(permutation_gml),
            "--verbosity",
            "fatal",
        ]
    )
    run_command(
        java_prefix
        + [
            "--input",
            str(input_path),
            "--matrixtype",
            "count",
            "--method",
            "ensemble",
            "--ensemblemethods",
            METHODS,
            "--ensembleparamfile",
            str(threshold_path),
            "--networkmergestrategy",
            "union",
            "--multigraph",
            "--filter",
            "rand",
            "--randroutine",
            "edgeScores",
            "--resamplemethod",
            "bootstrap",
            "--iterations",
            str(iterations),
            "--edgethreshold",
            "0.05",
            "--pvaluemerge",
            "brown",
            "--multicorr",
            "benjaminihochberg",
            "--nulldistribfile",
            str(permutation_path),
            "--randscorefile",
            str(bootstrap_path),
            "--scoreexport",
            "--format",
            "gml/tab_table",
            "--output",
            str(output_prefix),
            "--verbosity",
            "fatal",
        ]
    )
    return time.perf_counter() - start


def parse_conet_edges(tab_path: Path) -> pd.DataFrame:
    if not tab_path.exists() or tab_path.stat().st_size == 0:
        return pd.DataFrame(
            columns=["taxon_1", "taxon_2", "conet_selected", "conet_score", "conet_qvalue", "method_number"]
        )
    conet = pd.read_csv(tab_path, sep="\t", comment="#")
    if conet.empty:
        return pd.DataFrame(
            columns=["taxon_1", "taxon_2", "conet_selected", "conet_score", "conet_qvalue", "method_number"]
        )
    endpoints = conet["Label"].str.split("->", n=1, expand=True)
    conet[["taxon_1", "taxon_2"]] = endpoints.apply(
        lambda row: canonical_pair(row.iloc[0], row.iloc[1]),
        axis=1,
        result_type="expand",
    )
    conet = conet.drop_duplicates(["taxon_1", "taxon_2"]).copy()
    conet["conet_selected"] = 1
    qvalue = conet["qval"].clip(lower=1e-300)
    conet["conet_score"] = conet["method_number"].astype(float) + (-qvalue.map(math.log10) / 1000)
    conet["conet_qvalue"] = conet["qval"]
    return conet[["taxon_1", "taxon_2", "conet_selected", "conet_score", "conet_qvalue", "method_number"]]


def evaluate_dataset(
    dataset: str,
    n_samples: int,
    n_taxa: int,
    runtime_seconds: float,
    min_occurrence: int,
    iterations: int,
) -> pd.DataFrame:
    abundance = pd.read_csv(ROOT / "cleaned_data" / f"{dataset}_abundance.csv")
    eligible = set(abundance.columns[1:][(abundance.iloc[:, 1:] > 0).sum(axis=0) >= min_occurrence])

    truth = pd.read_csv(ROOT / "cleaned_data" / f"{dataset}_tested_pairs.csv")
    truth = truth[truth["tested_status"].astype(str).str.lower() != "ambiguous"].copy()
    truth[["taxon_1", "taxon_2"]] = truth.apply(
        lambda row: canonical_pair(row["taxon_1"], row["taxon_2"]),
        axis=1,
        result_type="expand",
    )
    truth = truth[truth["taxon_1"].isin(eligible) & truth["taxon_2"].isin(eligible)].copy()
    truth["interaction_label"] = truth["interaction_label"].astype(int)

    cache_dir = ROOT / "analysis_cache" / dataset / "conet"
    tab_path = cache_dir / f"conet_minocc{min_occurrence}_{iterations}.tab"
    conet = parse_conet_edges(tab_path)
    data = truth.merge(conet, on=["taxon_1", "taxon_2"], how="left", validate="one_to_one")
    data[["conet_selected", "conet_score", "method_number"]] = data[
        ["conet_selected", "conet_score", "method_number"]
    ].fillna(0)

    score_path = ROOT / "analysis_data" / f"{dataset}_conet_tested_pair_scores.csv"
    data[
        [
            "dataset",
            "taxon_1",
            "taxon_2",
            "interaction_label",
            "tested_status",
            "conet_selected",
            "conet_score",
            "conet_qvalue",
            "method_number",
        ]
    ].to_csv(score_path, index=False)

    labels = data["interaction_label"].astype(int)
    scores = data["conet_score"].astype(float)
    selected = data["conet_selected"].astype(int)
    positives = int(labels.sum())
    neutrals = int((labels == 0).sum())
    selected_edges = int(selected.sum())
    true_selected = int(((selected == 1) & (labels == 1)).sum())
    false_selected = int(((selected == 1) & (labels == 0)).sum())
    false_unselected = int(((selected == 0) & (labels == 1)).sum())

    precision = true_selected / (true_selected + false_selected) if true_selected + false_selected else math.nan
    recall = true_selected / (true_selected + false_unselected) if true_selected + false_unselected else math.nan
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall and not math.isnan(precision + recall)
        else math.nan
    )
    auroc = roc_auc_score(labels, scores) if positives and neutrals else math.nan
    average_precision = average_precision_score(labels, scores) if positives and neutrals else math.nan

    return pd.DataFrame(
        [
            {
                "dataset": dataset,
                "method": "CoNet",
                "min_occurrence": min_occurrence,
                "iterations": iterations,
                "n_samples": n_samples,
                "n_taxa": n_taxa,
                "tested_pairs": len(data),
                "positive_pairs": positives,
                "neutral_pairs": neutrals,
                "random_auprc": positives / len(data) if len(data) else math.nan,
                "average_precision": average_precision,
                "auroc": auroc,
                "native_selected_edges": selected_edges,
                "native_precision": precision,
                "native_recall": recall,
                "native_f1": f1,
                "runtime_seconds": runtime_seconds,
            }
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset")
    parser.add_argument("--min-occurrence", type=int, default=20)
    parser.add_argument("--iterations", type=int, default=100)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    cache_dir = ROOT / "analysis_cache" / args.dataset / "conet"
    cache_dir.mkdir(parents=True, exist_ok=True)
    n_samples, n_taxa, possible_pairs = prepare_input(args.dataset, cache_dir, args.min_occurrence)
    if n_taxa < 2:
        raise RuntimeError(f"{args.dataset}: fewer than two taxa survive the occurrence filter")

    # CoNet's top/bottom threshold guessing fails on small graphs if the
    # requested edge count makes the lower and upper thresholds overlap.
    max_edges = max(1, min(500, max(1, possible_pairs // 2)))
    runtime = run_conet(
        args.dataset,
        cache_dir,
        args.min_occurrence,
        args.iterations,
        max_edges,
        args.force,
    )
    result = evaluate_dataset(
        args.dataset,
        n_samples,
        n_taxa,
        runtime,
        args.min_occurrence,
        args.iterations,
    )
    output_path = ROOT / "results" / f"{args.dataset}_conet_performance.csv"
    output_path.parent.mkdir(exist_ok=True)
    result.to_csv(output_path, index=False)
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
