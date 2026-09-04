#!/usr/bin/env python3
"""Combine run counts and prepare the intact-community abundance table."""

import csv
from pathlib import Path
import re

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "external_data" / "carlstrom_phyllosphere_2019"
SOURCE = BASE / "source"
PROCESSED = BASE / "processed"
DATASET = "carlstrom_phyllosphere_2019"


def read_run_count(path):
    table = pd.read_csv(path, sep="\t", index_col=0)
    if table.shape[0] != 1:
        raise ValueError(f"Unexpected count table shape in {path}: {table.shape}")
    return table.iloc[0]


def main():
    report = pd.read_csv(BASE / "ena_read_runs.tsv", sep="\t")
    required_titles = set()
    metadata_tables = {}
    for experiment in (53, 61, 62):
        metadata = pd.read_csv(SOURCE / f"metadata{experiment}.csv", sep=";")
        metadata_tables[experiment] = metadata
        if experiment == 53:
            selected = metadata[
                (metadata.Treatment == "ALL") & metadata.Spray.isin(["U", "Mg"]) & (metadata.Time != "t0")
            ]
        else:
            selected = metadata[
                (~metadata.Treatment.isin(["Ax", "Axenic", "CORE"])) & (metadata.Time != "t0")
            ]
        required_titles.update(selected.Name)
    report = report[report.sample_title.isin(required_titles)]
    bulk_table = PROCESSED / "all_run_counts.tsv"
    if bulk_table.exists():
        raw = pd.read_csv(bulk_table, sep="\t", index_col=0)
        matrix = raw.T.fillna(0).astype(int)
        accession_to_title = dict(zip(report.run_accession, report.sample_title))
        matrix.index = [accession_to_title.get(value, value) for value in matrix.index]
        matrix.index.name = "sample_title"
    else:
        counts = []
        for row in report.itertuples(index=False):
            values = read_run_count(PROCESSED / "run_counts" / f"{row.run_accession}.tsv")
            values.name = row.sample_title
            counts.append(values)
        matrix = pd.DataFrame(counts).fillna(0).astype(int)
        matrix.index.name = "sample_title"
    matrix.to_csv(PROCESSED / "all_run_counts.csv")

    meta53 = metadata_tables[53]
    controls = meta53[
        (meta53.Treatment == "ALL") & meta53.Spray.isin(["U", "Mg"]) & (meta53.Time != "t0")
    ].copy()
    missing = sorted(set(controls.Name) - set(matrix.index))
    if missing:
        raise ValueError(f"Missing experiment 53 samples: {missing[:10]}")
    abundance = matrix.loc[controls.Name].copy()
    abundance = abundance.loc[:, abundance.sum(axis=0) > 0]
    abundance.insert(0, "sample_id", abundance.index)
    abundance.reset_index(drop=True).to_csv(ROOT / "cleaned_data" / f"{DATASET}_abundance.csv", index=False)

    removal_rows = []
    for experiment in (61, 62):
        metadata = metadata_tables[experiment]
        metadata = metadata[(~metadata.Treatment.isin(["Ax", "Axenic", "CORE"])) & (metadata.Time != "t0")]
        missing = sorted(set(metadata.Name) - set(matrix.index))
        if missing:
            raise ValueError(f"Missing experiment {experiment} samples: {missing[:10]}")
        for record in metadata.itertuples(index=False):
            row = matrix.loc[record.Name].to_dict()
            row.update({
                "sample_id": record.Name,
                "experiment": str(experiment),
                "treatment": record.Treatment,
                "replicate": record.Replicate,
            })
            removal_rows.append(row)
    removal = pd.DataFrame(removal_rows)
    identifiers = ["sample_id", "experiment", "treatment", "replicate"]
    removal[identifiers + [c for c in matrix.columns if c in removal.columns]].to_csv(
        PROCESSED / "removal_counts.csv", index=False
    )


if __name__ == "__main__":
    main()
