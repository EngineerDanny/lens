#!/usr/bin/env python3
"""Create analysis-ready abundance and tested-pair tables for Wortel SynCom."""

from pathlib import Path
import re

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "external_data/wortel_syncom_2026/source/SynCom"
OUT = ROOT / "cleaned_data"
DATASET = "wortel_syncom_2026"
TAXA = ["AM", "BA", "BOB", "BOV", "BT", "EC", "FP", "LB", "LJ"]
CUTOFF = 0.1


def members(text):
    return re.findall(r"'([^']+)'", str(text))


def solve_ct(equation, ct):
    match = re.fullmatch(r"\s*Y\s*=\s*([+-]?\d+(?:\.\d+)?)\s*\*\s*X\s*\+\s*([+-]?\d+(?:\.\d+)?)\s*", equation)
    if not match:
        raise ValueError(f"Unsupported standard curve: {equation}")
    slope, intercept = map(float, match.groups())
    return max(0.0, (float(ct) - intercept) / slope)


def wide_replicates(frame, group_col, taxon_col, replicate_cols, prefix):
    rows = []
    for group, part in frame.groupby(group_col, sort=False):
        for replicate, column in enumerate(replicate_cols, 1):
            row = {taxon: 0.0 for taxon in TAXA}
            for _, record in part.iterrows():
                taxon = record[taxon_col]
                if taxon in row:
                    row[taxon] = solve_ct(record["equation"], record[column])
            row["sample_id"] = f"{prefix}_{group}_{replicate}"
            rows.append(row)
    return rows


def prepare_abundance():
    pair = pd.read_csv(SOURCE / "QPCR data/Files/In_vitro_data_with_raw.csv")
    pair["members"] = pair["culture"].map(members)
    pair = pair[pair["members"].map(lambda x: "RI" not in x)].copy()
    rows = wide_replicates(pair, "culture", "focal_species", ["Bio1", "Bio2", "Bio3", "Bio4"], "pair")

    higher = pd.read_csv(SOURCE / "Main/Input Files/Culture Experiments/in_vitro_df_triquad_with_solutions.csv")
    rows.extend(wide_replicates(higher, "Community", "focal_spec", ["Bio1_avg", "Bio2_avg"], "higher_a"))

    second = pd.read_csv(SOURCE / "Main/Input Files/Culture Experiments/in_vitro_df_second_with_solutions.csv")
    rows.extend(wide_replicates(second, "Community", "Species", ["Bio1", "Bio2"], "higher_b"))

    abundance = pd.DataFrame(rows)[["sample_id", *TAXA]]
    if abundance["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample identifiers")
    if not np.isfinite(abundance[TAXA].to_numpy()).all() or (abundance[TAXA] < 0).any().any():
        raise ValueError("Invalid abundance value")
    if (abundance[TAXA].sum(axis=1) <= 0).any():
        raise ValueError("Empty abundance sample")
    abundance.to_csv(OUT / f"{DATASET}_abundance.csv", index=False)
    return abundance


def classify(value):
    if value > CUTOFF:
        return 1
    if value < -CUTOFF:
        return -1
    return 0


def prepare_truth():
    effects = pd.read_csv(SOURCE / "Parameters/output/epsilons.csv")
    directed = {
        (row.bacteria, row.medium): {
            "strength": row.eps_yield,
            "class": classify(row.eps_yield),
        }
        for row in effects.itertuples()
    }
    no_growth = [("AM", "BA"), ("AM", "BT"), ("FP", "BA"), ("FP", "BOB"), ("FP", "BOV"), ("FP", "BT")]
    for key in no_growth:
        directed[key] = {"strength": np.nan, "class": -1}

    rows = []
    for index, taxon_1 in enumerate(TAXA):
        for taxon_2 in TAXA[index + 1:]:
            first = directed[(taxon_1, taxon_2)]
            second = directed[(taxon_2, taxon_1)]
            classes = (first["class"], second["class"])
            interaction = int(classes != (0, 0))
            nonzero = {value for value in classes if value != 0}
            effect_sign = next(iter(nonzero)) if len(nonzero) == 1 else np.nan
            strengths = [abs(x["strength"]) for x in (first, second) if np.isfinite(x["strength"])]
            rows.append({
                "dataset": DATASET,
                "taxon_1": taxon_1,
                "taxon_2": taxon_2,
                "tested_status": "positive" if interaction else "neutral",
                "interaction_label": interaction,
                "effect_sign": effect_sign,
                "effect_strength": max(strengths) if strengths else np.nan,
                "n_evidence": 2,
                "truth_type": "conditioned_medium_final_population_density",
                "label_rule": "edge if either directed final-population-density log ratio is outside [-0.1, 0.1]; no-growth outcomes are negative effects",
                "mixed_evidence": len(nonzero) > 1,
                "direction_available": True,
                "experimental_setting": "conditioned YCFA medium and separate qPCR cocultures",
            })
    truth = pd.DataFrame(rows)
    truth.to_csv(OUT / f"{DATASET}_tested_pairs.csv", index=False)
    return truth


def main():
    OUT.mkdir(exist_ok=True)
    abundance = prepare_abundance()
    truth = prepare_truth()
    manifest_path = OUT / "manifest.csv"
    manifest = pd.read_csv(manifest_path)
    manifest = manifest[manifest["dataset"] != DATASET]
    manifest = pd.concat([manifest, pd.DataFrame([{
        "dataset": DATASET,
        "n_samples": len(abundance),
        "n_taxa": len(TAXA),
        "zero_frequency": (abundance[TAXA] == 0).to_numpy().mean(),
        "n_tested_pairs": len(truth),
        "n_positive": int(truth["interaction_label"].sum()),
        "n_tested_neutral": int((truth["interaction_label"] == 0).sum()),
        "n_ambiguous": 0,
        "n_excluded_records": 9,
        "source_abundance": "external_data/wortel_syncom_2026/source/SynCom/QPCR data/Files/In_vitro_data_with_raw.csv",
    }])], ignore_index=True)
    manifest.to_csv(manifest_path, index=False)
    print(f"abundance: {abundance.shape[0]} samples x {abundance.shape[1] - 1} taxa")
    print(f"truth: {len(truth)} tested pairs; {truth.interaction_label.sum()} interactions; {(truth.interaction_label == 0).sum()} neutrals")
    print(f"zero frequency: {(abundance[TAXA] == 0).to_numpy().mean():.4f}")


if __name__ == "__main__":
    main()
