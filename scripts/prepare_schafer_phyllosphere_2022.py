#!/usr/bin/env python3
"""Prepare the Schäfer et al. phyllosphere perturbation benchmark."""

from pathlib import Path
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "external_data/schafer_phyllosphere_2022/source"
OUT = ROOT / "cleaned_data"
DATASET = "schafer_phyllosphere_2022"


def main():
    counts = pd.read_excel(SOURCE / "supplementary_data_1.xlsx", sheet_name="Data_table 3", header=2, index_col=0)
    metadata = pd.read_excel(SOURCE / "supplementary_data_1.xlsx", sheet_name="Data_table 4", header=2)
    metadata = metadata.rename(columns={metadata.columns[0]: "sample_id"}).set_index("sample_id")
    truth = pd.read_excel(SOURCE / "source_data_extended_figure_1.xlsx")

    counts = counts.drop(index="Unclassified", errors="ignore").apply(pd.to_numeric, errors="coerce").fillna(0)
    sample_ids = metadata.index[metadata["Treatment"].isin(["Screen", "FocalCom"])]
    abundance = counts.loc[:, counts.columns.intersection(sample_ids)].T
    abundance = abundance.loc[:, abundance.nunique() > 1]
    abundance.index.name = "sample_id"
    abundance.reset_index().to_csv(OUT / f"{DATASET}_abundance.csv", index=False)

    mapping = metadata[["Condition", "Added_ASV"]].dropna().drop_duplicates()
    if (mapping.groupby("Condition")["Added_ASV"].nunique() > 1).any():
        raise ValueError("An added strain maps to multiple ASVs")
    truth["added_asv"] = truth["Added_strain"].map(dict(zip(mapping["Condition"], mapping["Added_ASV"])))
    truth = truth[truth["significant"].notna() & truth["added_asv"].notna()].copy()
    truth = truth[truth["Focal_strain"].isin(abundance.columns) & truth["added_asv"].isin(abundance.columns)]
    truth["taxon_1"] = truth[["Focal_strain", "added_asv"]].min(axis=1)
    truth["taxon_2"] = truth[["Focal_strain", "added_asv"]].max(axis=1)
    truth = truth[truth["taxon_1"] != truth["taxon_2"]]

    rows = []
    for (taxon_1, taxon_2), group in truth.groupby(["taxon_1", "taxon_2"], sort=True):
        labels = set(group["significant"].astype(int))
        if len(labels) != 1:
            continue
        label = labels.pop()
        significant = group[group["significant"] == 1]
        signs = set(np.sign(significant["log2FC"]).astype(int)) if len(significant) else set()
        rows.append({
            "dataset": DATASET,
            "taxon_1": taxon_1,
            "taxon_2": taxon_2,
            "tested_status": "positive" if label else "neutral",
            "interaction_label": label,
            "effect_sign": next(iter(signs)) if len(signs) == 1 else np.nan,
            "effect_strength": significant["log2FC"].abs().max() if len(significant) else group["log2FC"].abs().max(),
            "n_evidence": len(group),
            "truth_type": "strain_addition_effect_on_focal_community",
            "label_rule": "edge when adjusted Wald-test p <= 0.01; pairs with conflicting strain evidence after ASV mapping are excluded",
            "mixed_evidence": False,
            "direction_available": True,
            "experimental_setting": "Arabidopsis phyllosphere 15-strain focal community with one added strain",
        })
    tested = pd.DataFrame(rows)
    tested.to_csv(OUT / f"{DATASET}_tested_pairs.csv", index=False)

    manifest_path = OUT / "manifest.csv"
    manifest = pd.read_csv(manifest_path)
    manifest = manifest[manifest["dataset"] != DATASET]
    manifest = pd.concat([manifest, pd.DataFrame([{
        "dataset": DATASET,
        "n_samples": len(abundance),
        "n_taxa": abundance.shape[1],
        "zero_frequency": (abundance.to_numpy() == 0).mean(),
        "n_tested_pairs": len(tested),
        "n_positive": int(tested["interaction_label"].sum()),
        "n_tested_neutral": int((tested["interaction_label"] == 0).sum()),
        "n_ambiguous": 25,
        "n_excluded_records": 80,
        "source_abundance": "external_data/schafer_phyllosphere_2022/source/supplementary_data_1.xlsx",
    }])], ignore_index=True)
    manifest.to_csv(manifest_path, index=False)

    print(f"abundance: {len(abundance)} samples x {abundance.shape[1]} taxa")
    print(f"truth: {len(tested)} tested pairs; {tested.interaction_label.sum()} interactions; {(tested.interaction_label == 0).sum()} neutrals")
    print(f"zero frequency: {(abundance.to_numpy() == 0).mean():.4f}")


if __name__ == "__main__":
    main()
