#!/usr/bin/env python3
"""Convert a cleaned sample-by-taxon CSV to CoNet's taxa-by-sample TSV."""

import argparse
import csv


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_csv")
    parser.add_argument("output_tsv")
    parser.add_argument("--min-occurrence", type=int, default=0)
    args = parser.parse_args()

    with open(args.input_csv, newline="") as handle:
        rows = list(csv.reader(handle))
    header = rows[0]
    sample_ids = [row[0] for row in rows[1:]]

    with open(args.output_tsv, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon_id", *sample_ids])
        for column, taxon in enumerate(header[1:], start=1):
            values = [row[column] for row in rows[1:]]
            if sum(float(value) > 0 for value in values) >= args.min_occurrence:
                writer.writerow([taxon, *values])


if __name__ == "__main__":
    main()
