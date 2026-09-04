#!/usr/bin/env python3
"""Join pair features to labels and create fixed study holdout folds."""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CLEAN = ROOT / "cleaned_data"
FEATURES = ROOT / "analysis_data"
OUTPUT = ROOT / "training_data"

DATASETS = [
    "omm12",
    "omm12_keystone_2023",
    "pairinterax",
    "butyrate_assembly_2021",
    "host_fitness_2018",
]

FOLDS = [
    ("fold1", ["omm12"]),
    ("fold2", ["omm12_keystone_2023"]),
    ("fold3", ["pairinterax"]),
    ("fold4", ["butyrate_assembly_2021"]),
    ("fold5", ["host_fitness_2018"]),
    ("fold6", ["omm12", "omm12_keystone_2023"]),
]

LABEL_FIELDS = [
    "tested_status",
    "interaction_label",
    "effect_sign",
    "effect_strength",
    "n_evidence",
    "truth_type",
    "label_rule",
    "mixed_evidence",
    "direction_available",
    "experimental_setting",
]

PRIMARY_PREDICTORS = [
    "pln_score_percentile",
    "pln_selected",
    "glmnet_score_percentile",
    "glmnet_selected",
    "spieceasi_score_percentile",
    "spieceasi_selected",
    "spieceasi_stability",
    "spring_score_percentile",
    "spring_selected",
    "spring_stability",
    "sparcc_score_percentile",
    "sparcc_selected",
    "prevalence_1",
    "prevalence_2",
    "prevalence_min",
    "prevalence_max",
    "joint_prevalence",
    "n_joint_positive",
    "presence_jaccard",
    "n_samples",
    "n_taxa",
    "overall_zero_frequency",
]


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def key(row: dict[str, str]) -> tuple[str, str, str]:
    taxon_1, taxon_2 = sorted((row["taxon_1"], row["taxon_2"]))
    return row["dataset"], taxon_1, taxon_2


def main() -> None:
    joined: list[dict[str, str]] = []
    feature_fields: list[str] | None = None

    for dataset in DATASETS:
        fields, feature_rows = read_csv(FEATURES / f"{dataset}_pair_features.csv")
        current_feature_fields = [field for field in fields if field not in {"dataset", "taxon_1", "taxon_2"}]
        if feature_fields is None:
            feature_fields = current_feature_fields
        elif current_feature_fields != feature_fields:
            raise ValueError(f"Feature schema differs for {dataset}")

        _, labels = read_csv(CLEAN / f"{dataset}_tested_pairs.csv")
        feature_map = {key(row): row for row in feature_rows}
        if len(feature_map) != len(feature_rows):
            raise ValueError(f"Duplicate feature keys in {dataset}")

        for label in labels:
            pair = key(label)
            if pair not in feature_map:
                raise ValueError(f"Missing feature row for {pair}")
            row = {
                "dataset": dataset,
                "taxon_1": label["taxon_1"],
                "taxon_2": label["taxon_2"],
            }
            row.update({field: label.get(field, "") for field in LABEL_FIELDS})
            row.update({field: feature_map[pair].get(field, "") for field in feature_fields})
            joined.append(row)

    assert feature_fields is not None
    fields = ["dataset", "taxon_1", "taxon_2", *LABEL_FIELDS, *feature_fields]
    joined.sort(key=lambda row: (DATASETS.index(row["dataset"]), row["taxon_1"], row["taxon_2"]))

    write_csv(OUTPUT / "all_tested_pairs.csv", fields, joined)
    primary = [row for row in joined if row["tested_status"] in {"positive", "tested_neutral"}]
    for row in primary:
        for predictor in PRIMARY_PREDICTORS:
            value = row.get(predictor, "")
            if value == "" or not math.isfinite(float(value)):
                raise ValueError(f"Invalid primary predictor {predictor} for {key(row)}")
    write_csv(OUTPUT / "all_labeled_pairs.csv", fields, primary)
    write_csv(
        OUTPUT / "model_columns.csv",
        ["column", "role"],
        [{"column": column, "role": "primary_predictor"} for column in PRIMARY_PREDICTORS]
        + [{"column": "interaction_label", "role": "response"}],
    )

    manifest_rows = []
    for fold_name, held_out in FOLDS:
        train = [row for row in primary if row["dataset"] not in held_out]
        test = [row for row in primary if row["dataset"] in held_out]
        if set(row["dataset"] for row in train) & set(row["dataset"] for row in test):
            raise AssertionError(f"Study leakage in {fold_name}")
        fold_dir = OUTPUT / fold_name
        write_csv(fold_dir / "train.csv", fields, train)
        write_csv(fold_dir / "test.csv", fields, test)
        manifest_rows.append({
            "fold": fold_name,
            "held_out_studies": ";".join(held_out),
            "n_train": len(train),
            "n_train_positive": sum(row["interaction_label"] == "1" for row in train),
            "n_train_neutral": sum(row["interaction_label"] == "0" for row in train),
            "n_test": len(test),
            "n_test_positive": sum(row["interaction_label"] == "1" for row in test),
            "n_test_neutral": sum(row["interaction_label"] == "0" for row in test),
        })

    write_csv(
        OUTPUT / "fold_manifest.csv",
        [
            "fold", "held_out_studies", "n_train", "n_train_positive",
            "n_train_neutral", "n_test", "n_test_positive", "n_test_neutral",
        ],
        manifest_rows,
    )


if __name__ == "__main__":
    main()
