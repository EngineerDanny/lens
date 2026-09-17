#!/usr/bin/env python3
"""Create analysis-ready abundance and tested-pair tables.

The script reads the source files under interaction_ground_truth and writes a
separate cleaned_data tree. Source files are never modified.
"""

from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "interaction_ground_truth"
OUTPUT = ROOT / "cleaned_data"

DATASETS = {"butyrate_assembly_2021": "abundance_matrix.tsv.gz"}


PAIR_FIELDS = [
    "dataset",
    "taxon_1",
    "taxon_2",
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


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_csv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def numeric(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def pair_key(a: str, b: str) -> tuple[str, str]:
    return tuple(sorted((a.strip(), b.strip())))


def bool_text(value: bool) -> str:
    return "TRUE" if value else "FALSE"


def consensus_sign(values: list[int]) -> str:
    nonzero = {value for value in values if value != 0}
    return str(next(iter(nonzero))) if len(nonzero) == 1 else ""


def clean_abundance(dataset: str) -> tuple[list[str], list[dict[str, str]], list[dict[str, str]]]:
    source_path = SOURCE / dataset / "processed" / DATASETS[dataset]
    fields, rows = read_tsv_gz(source_path)
    exclusions: list[dict[str, str]] = []
    return fields, rows, exclusions


def threshold_pairs(dataset: str, truth_type: str, context: str) -> list[dict[str, object]]:
    proc = SOURCE / dataset / "processed"
    _, raw = read_tsv_gz(proc / "truth_edges_raw.tsv.gz")
    _, positives = read_tsv_gz(proc / "truth_undirected.tsv.gz")
    positive_map = {
        pair_key(row["taxon_1"], row["taxon_2"]): row
        for row in positives
    }
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        a, b = row["partner_taxon"], row["focal_taxon"]
        if a != b:
            grouped[pair_key(a, b)].append(row)

    result = []
    for (a, b), evidence in sorted(grouped.items()):
        positive = (a, b) in positive_map
        truth = positive_map.get((a, b), {})
        signs = [int(row["effect_sign"]) for row in evidence if row["effect_sign"]]
        contexts = sorted({row.get("context", "") for row in evidence if row.get("context", "")})
        result.append({
            "dataset": dataset,
            "taxon_1": a,
            "taxon_2": b,
            "tested_status": "positive" if positive else "tested_neutral",
            "interaction_label": 1 if positive else 0,
            "effect_sign": truth.get("sign_consensus", "") if positive else consensus_sign(signs),
            "effect_strength": truth.get("max_abs_strength", "") if positive else "",
            "n_evidence": len(evidence),
            "truth_type": truth_type,
            "label_rule": "positive when absolute aggregated log2 ratio is at least 0.5",
            "mixed_evidence": bool_text(len(set(signs)) > 1),
            "direction_available": "TRUE",
            "experimental_setting": ";".join(contexts) or context,
        })
    return result


def validate_abundance(fields: list[str], rows: list[dict[str, str]]) -> dict[str, object]:
    taxa = fields[1:]
    assert fields[0] == "sample_id"
    assert len({row["sample_id"] for row in rows}) == len(rows)
    assert len(set(taxa)) == len(taxa)
    values = []
    for row in rows:
        for taxon in taxa:
            value = numeric(row[taxon])
            assert value is not None and value >= 0
            values.append(value)
    assert all(any(numeric(row[taxon]) != 0 for row in rows) for taxon in taxa)
    return {
        "n_samples": len(rows),
        "n_taxa": len(taxa),
        "zero_frequency": sum(value == 0 for value in values) / len(values),
    }


def validate_pairs(rows: list[dict[str, object]], taxa: set[str]) -> dict[str, int]:
    keys = [(str(row["taxon_1"]), str(row["taxon_2"])) for row in rows]
    assert len(keys) == len(set(keys))
    assert all(a < b for a, b in keys)
    assert all(a in taxa and b in taxa for a, b in keys)
    allowed = {"positive", "tested_neutral", "ambiguous"}
    assert all(row["tested_status"] in allowed for row in rows)
    assert all(
        row["interaction_label"] in (1, "1") if row["tested_status"] == "positive"
        else row["interaction_label"] in (0, "0") if row["tested_status"] == "tested_neutral"
        else row["interaction_label"] == ""
        for row in rows
    )
    counts = defaultdict(int)
    for row in rows:
        counts[str(row["tested_status"])] += 1
    return dict(counts)


def main() -> None:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    all_exclusions: list[dict[str, str]] = []
    manifest_path = OUTPUT / "manifest.csv"
    manifest = []
    if manifest_path.exists():
        with manifest_path.open(newline="") as handle:
            manifest = [row for row in csv.DictReader(handle) if row["dataset"] not in DATASETS]

    for dataset in DATASETS:
        fields, abundance, exclusions = clean_abundance(dataset)
        taxa = fields[1:]
        taxa_set = set(taxa)

        pairs = threshold_pairs(dataset, "community_assembly", "defined-community assembly")

        profile = validate_abundance(fields, abundance)
        label_counts = validate_pairs(pairs, taxa_set)

        write_csv(OUTPUT / f"{dataset}_abundance.csv", fields, abundance)
        write_csv(OUTPUT / f"{dataset}_tested_pairs.csv", PAIR_FIELDS, pairs)

        all_exclusions.extend(exclusions)
        manifest.append({
            "dataset": dataset,
            **profile,
            "n_tested_pairs": len(pairs),
            "n_positive": label_counts.get("positive", 0),
            "n_tested_neutral": label_counts.get("tested_neutral", 0),
            "n_ambiguous": label_counts.get("ambiguous", 0),
            "source_abundance": str((SOURCE / dataset / "processed" / DATASETS[dataset]).relative_to(ROOT)),
        })

    write_csv(
        OUTPUT / "manifest.csv",
        [
            "dataset", "n_samples", "n_taxa", "zero_frequency", "n_tested_pairs",
            "n_positive", "n_tested_neutral", "n_ambiguous", "n_excluded_records",
            "source_abundance",
        ],
        [
            {
                **row,
                "n_excluded_records": row.get("n_excluded_records", sum(item["dataset"] == row["dataset"] for item in all_exclusions)),
            }
            for row in manifest
        ],
    )


if __name__ == "__main__":
    main()
