#!/usr/bin/env python3
"""Prepare the Friedman, Higgins and Gore 2017 microcosm benchmark."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
DATASET = "friedman_microcosm_2017"
SOURCE = (
    ROOT
    / "external_data"
    / DATASET
    / "source"
    / "jfriedman-micro-micro_assembly_rules-7adbaa1"
    / "data"
)
TAXA = ["Ea", "Pa", "Pch", "Pci", "Pf", "Pp", "Pv", "Sm"]


def endpoint_rows(path: Path) -> list[dict]:
    records: list[dict] = []
    workbook = pd.ExcelFile(path)
    for sheet in workbook.sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet)
        time = pd.to_numeric(frame.iloc[:, 0], errors="coerce")
        endpoints = frame.loc[time == 5, frame.columns[1:]]
        for replicate, (_, row) in enumerate(endpoints.iterrows(), start=1):
            record = {taxon: 0.0 for taxon in TAXA}
            for taxon in endpoints.columns:
                record[taxon] = float(row[taxon]) if pd.notna(row[taxon]) else np.nan
            if any(pd.isna(record[taxon]) for taxon in endpoints.columns):
                continue
            record["sample_id"] = f"{path.stem}__{sheet}__endpoint_{replicate:02d}"
            records.append(record)
    return records


def raw_endpoint_count(path: Path) -> int:
    count = 0
    for sheet in pd.ExcelFile(path).sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet)
        count += int((pd.to_numeric(frame.iloc[:, 0], errors="coerce") == 5).sum())
    return count


def monoculture_endpoints() -> dict[str, np.ndarray]:
    path = SOURCE / "monoculture_timeSeries.xlsx"
    result = {}
    for sheet in pd.ExcelFile(path).sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet)
        time = pd.to_numeric(frame.iloc[:, 0], errors="coerce")
        values = pd.to_numeric(frame.loc[time == 5].iloc[0, 1:], errors="coerce").dropna()
        result[sheet] = values.to_numpy(float)
    return result


def directional_truth() -> pd.DataFrame:
    monoculture = monoculture_endpoints()
    path = SOURCE / "pair_timeSeries.xlsx"
    rows = []
    for sheet in pd.ExcelFile(path).sheet_names:
        frame = pd.read_excel(path, sheet_name=sheet)
        time = pd.to_numeric(frame.iloc[:, 0], errors="coerce")
        endpoints = frame.loc[time == 5, frame.columns[1:]]
        pair_taxa = list(endpoints.columns)
        for focal in pair_taxa:
            partner = pair_taxa[1] if focal == pair_taxa[0] else pair_taxa[0]
            pair_values = pd.to_numeric(endpoints[focal], errors="coerce").dropna().to_numpy(float)
            mono_values = monoculture[focal]
            _, pvalue = ttest_ind(pair_values, mono_values, equal_var=False)
            pair_mean = float(pair_values.mean())
            mono_mean = float(mono_values.mean())
            relative_yield = (pair_mean - mono_mean) / (pair_mean + mono_mean)
            rows.append({
                "dataset": DATASET,
                "experimental_study": "Friedman_Higgins_Gore_2017",
                "taxon_1": focal,
                "taxon_2": partner,
                "direction": f"{partner}_to_{focal}",
                "effect_strength": relative_yield,
                "effect_sign": int(np.sign(relative_yield)),
                "n_pair_replicates": len(pair_values),
                "n_monoculture_replicates": len(mono_values),
                "p_value": pvalue,
            })
    truth = pd.DataFrame(rows)
    truth["q_value"] = multipletests(truth["p_value"], method="fdr_bh")[1]
    truth["interaction_label"] = (truth["q_value"] <= 0.05).astype(int)
    truth["tested_status"] = np.where(truth["interaction_label"] == 1, "positive", "neutral")
    truth["truth_type"] = "directional_pairwise_coculture_effect"
    truth["label_rule"] = "Welch test of endpoint abundance against monoculture; BH q <= 0.05"
    truth["direction_available"] = True
    truth["experimental_setting"] = "soil bacterial microcosm; five growth cycles"
    return truth.sort_values(["taxon_1", "taxon_2"]).reset_index(drop=True)


def undirected_truth(directional: pd.DataFrame) -> pd.DataFrame:
    frame = directional.copy()
    ordered = np.sort(frame[["taxon_1", "taxon_2"]].astype(str).to_numpy(), axis=1)
    frame[["taxon_1", "taxon_2"]] = ordered
    grouped = frame.groupby(["dataset", "taxon_1", "taxon_2"], as_index=False)
    output = grouped.agg(
        interaction_label=("interaction_label", "max"),
        effect_strength=("effect_strength", lambda x: float(x.iloc[np.argmax(np.abs(x.to_numpy()))])),
        n_evidence=("interaction_label", "size"),
    )
    output["tested_status"] = np.where(output["interaction_label"] == 1, "positive", "neutral")
    output["effect_sign"] = np.sign(output["effect_strength"]).astype(int)
    output["truth_type"] = "undirected_summary_of_directional_pairwise_effects"
    output["label_rule"] = "positive when either tested direction has BH q <= 0.05"
    output["direction_available"] = True
    output["experimental_setting"] = "soil bacterial microcosm; five growth cycles"
    return output


def main() -> None:
    abundance_paths = [
        SOURCE / "trio_lastTransfer.xlsx",
        SOURCE / "7and8Species_lastTransfer.xlsx",
    ]
    raw_endpoints = sum(raw_endpoint_count(path) for path in abundance_paths)
    abundance_rows = endpoint_rows(abundance_paths[0])
    abundance_rows += endpoint_rows(abundance_paths[1])
    abundance = pd.DataFrame(abundance_rows)[["sample_id"] + TAXA]
    incomplete_endpoints = raw_endpoints - len(abundance)
    zero_total = abundance[TAXA].sum(axis=1) <= 0
    zero_total_endpoints = int(zero_total.sum())
    abundance = abundance.loc[~zero_total].reset_index(drop=True)
    if abundance.isna().any().any() or (abundance[TAXA].sum(axis=1) <= 0).any():
        raise RuntimeError("Invalid higher-order abundance table")

    directional = directional_truth()
    undirected = undirected_truth(directional)
    cleaned = ROOT / "cleaned_data"
    cleaned.mkdir(exist_ok=True)
    abundance.to_csv(cleaned / f"{DATASET}_abundance.csv", index=False)
    directional.to_csv(cleaned / f"{DATASET}_tested_pairs_directional.csv", index=False)
    undirected.to_csv(cleaned / f"{DATASET}_tested_pairs.csv", index=False)

    audit = {
        "dataset": DATASET,
        "abundance_source": "trio and seven/eight-species endpoint cultures only",
        "truth_source": "monoculture and pair-culture endpoint experiments only",
        "n_samples": len(abundance),
        "excluded_incomplete_endpoint_samples": incomplete_endpoints,
        "excluded_zero_total_endpoint_samples": zero_total_endpoints,
        "n_taxa": len(TAXA),
        "directed_tested_effects": len(directional),
        "directed_supported_effects": int(directional["interaction_label"].sum()),
        "directed_tested_neutrals": int((directional["interaction_label"] == 0).sum()),
        "undirected_tested_pairs": len(undirected),
        "undirected_positive_pairs": int(undirected["interaction_label"].sum()),
        "undirected_neutral_pairs": int((undirected["interaction_label"] == 0).sum()),
        "network_truth_experiment_overlap": False,
    }
    output = ROOT / "results" / f"{DATASET}_audit"
    output.mkdir(parents=True, exist_ok=True)
    with open(output / "audit_summary.json", "w") as handle:
        json.dump(audit, handle, indent=2)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
