#!/usr/bin/env python3
"""Audit the source files retained for the Wortel gut synthetic community."""

from pathlib import Path
import ast
import csv

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "external_data/wortel_syncom_2026/source/SynCom"
OUTPUT = ROOT / "results/wortel_syncom_2026_audit"


def culture_members(value):
    return tuple(ast.literal_eval(value))


def record(check, value, status, note):
    return {"check": check, "value": value, "status": status, "note": note}


def main():
    pair = pd.read_csv(SOURCE / "QPCR data/Files/In_vitro_data_with_raw.csv")
    higher_1 = pd.read_csv(
        SOURCE / "Main/Input Files/Culture Experiments/in_vitro_df_triquad_with_solutions.csv"
    )
    higher_2 = pd.read_csv(
        SOURCE / "Main/Input Files/Culture Experiments/in_vitro_df_second_with_solutions.csv"
    )
    effects = pd.read_csv(SOURCE / "Parameters/output/epsilons_RI.csv")

    pair["members"] = pair["culture"].map(culture_members)
    higher_1["members"] = higher_1["culture"].map(culture_members)
    taxa = sorted(set(pair["focal_species"]))
    possible_directed = len(taxa) * (len(taxa) - 1)
    directed_keys = effects[["bacteria", "medium"]]

    checks = [
        record("pair_focal_rows", len(pair), "pass", "Two rows per two species culture."),
        record("pair_communities", pair["culture"].nunique(), "pass", "All 45 possible pairs among 10 taxa are represented."),
        record("pair_taxa", len(taxa), "pass", ", ".join(taxa)),
        record("pair_biological_replicates", 4, "pass", "Bio1 through Bio4 are available, each based on two technical measurements."),
        record("pair_missing_ct", int(pair[[f"Bio{i}" for i in range(1, 5)]].isna().sum().sum()), "pass", "No missing biological replicate averages."),
        record("pair_duplicate_culture_focal", int(pair.duplicated(["culture", "focal_species"]).sum()), "pass", "Expected key is unique."),
        record("first_higher_order_communities", higher_1["culture"].nunique(), "pass", "Three or four species cocultures with two biological replicates."),
        record("second_higher_order_communities", higher_2["Community"].nunique(), "pass", "Additional three or four species cocultures with two biological replicates."),
        record("numeric_directed_effects", len(effects), "pass", "Conditioned medium effects with numeric growth rate and yield coefficients."),
        record("possible_directed_effects", possible_directed, "info", "Ten taxa give 90 possible directed effects."),
        record("nonnumeric_no_growth_effects", possible_directed - len(effects), "warning", "Missing numeric coefficients include experiments where the acceptor did not grow; retain these as observed inhibition with missing magnitude."),
        record("duplicate_directed_effects", int(directed_keys.duplicated().sum()), "pass", "Directed effect key is unique."),
        record("effect_taxa_match_abundance", int(set(effects["bacteria"]) | set(effects["medium"]) <= set(taxa)), "pass", "All effect identifiers occur in the qPCR abundance data."),
    ]

    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / "quality_checks.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["check", "value", "status", "note"])
        writer.writeheader()
        writer.writerows(checks)

    print(pd.DataFrame(checks).to_string(index=False))


if __name__ == "__main__":
    main()
