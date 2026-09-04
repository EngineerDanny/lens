#!/usr/bin/env python3
"""Process the required Carlström reads in one multithreaded VSEARCH run."""

import csv
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "external_data" / "carlstrom_phyllosphere_2019"
READS = BASE / "raw_reads"
SOURCE = BASE / "source"
OUT = BASE / "processed"
REPORT = BASE / "ena_read_runs.tsv"
REFERENCE = SOURCE / "64_strains_sequences.fasta"


def required_titles():
    keep = set()
    for experiment in (53, 61, 62):
        with (SOURCE / f"metadata{experiment}.csv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter=";"))
        for row in rows:
            if experiment == 53:
                selected = row["Treatment"] == "ALL" and row["Spray"] in {"U", "Mg"} and row["Time"] != "t0"
            else:
                selected = row["Treatment"] not in {"Ax", "Axenic", "CORE"} and row["Time"] != "t0"
            if selected:
                keep.add(row["Name"])
    return keep


def combine(rows, mate, destination):
    with destination.open("wb") as output:
        for row in rows:
            path = READS / f"{row['run_accession']}_{mate}.fastq.gz"
            with path.open("rb") as handle:
                shutil.copyfileobj(handle, output, length=1024 * 1024)


def run(command):
    subprocess.run(command, check=True)


def split_fasta(filtered, chunk_paths):
    handles = [path.open("w") for path in chunk_paths]
    try:
        index = -1
        output = None
        with filtered.open() as source:
            for line in source:
                if line.startswith(">"):
                    index += 1
                    output = handles[index % len(handles)]
                output.write(line)
    finally:
        for handle in handles:
            handle.close()


def map_chunk(chunk):
    counts = Counter()
    mapper = subprocess.Popen(
        ["minimap2", "-x", "sr", "-t", "1", "--secondary=no", str(REFERENCE), str(chunk)],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
    )
    awk_program = (
        'BEGIN { FS="\\t"; OFS="\\t" } '
        '$11 > 0 && $10 / $11 >= 0.97 && $11 / $2 >= 0.90 '
        '{ split($1, a, "."); key=$6 SUBSEP a[1]; count[key]++ } '
        'END { for (key in count) { split(key, value, SUBSEP); print value[1], value[2], count[key] } }'
    )
    aggregator = subprocess.Popen(
        ["awk", awk_program], stdin=mapper.stdout, stdout=subprocess.PIPE, text=True,
    )
    assert mapper.stdout is not None
    mapper.stdout.close()
    output, _ = aggregator.communicate()
    if mapper.wait() != 0 or aggregator.returncode != 0:
        raise RuntimeError("minimap2 mapping failed")
    for line in output.splitlines():
        target, accession, value = line.split("\t")
        counts[(target, accession)] += int(value)
    return counts


def map_with_minimap2(filtered, table_path, rows):
    chunk_paths = [OUT / f"filtered_chunk_{index}.fasta" for index in range(8)]
    if not all(path.exists() for path in chunk_paths):
        for path in chunk_paths:
            path.unlink(missing_ok=True)
        split_fasta(filtered, chunk_paths)
    counts = Counter()
    with ThreadPoolExecutor(max_workers=8) as executor:
        for result in executor.map(map_chunk, chunk_paths):
            counts.update(result)
    accessions = [row["run_accession"] for row in rows]
    targets = []
    with REFERENCE.open() as handle:
        for line in handle:
            if line.startswith(">"):
                targets.append(line[1:].split()[0])
    with table_path.open("w") as output:
        output.write("#OTU ID\t" + "\t".join(accessions) + "\n")
        for target in targets:
            values = [str(counts[(target, accession)]) for accession in accessions]
            output.write(target + "\t" + "\t".join(values) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    with REPORT.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["sample_title"] in required_titles()]
    if len(rows) != 670:
        raise ValueError(f"Expected 670 required runs, found {len(rows)}")
    forward = OUT / "required_R1.fastq.gz"
    reverse = OUT / "required_R2.fastq.gz"
    if not forward.exists():
        combine(rows, 1, forward)
    if not reverse.exists():
        combine(rows, 2, reverse)
    merged = OUT / "required_merged.fastq"
    filtered = OUT / "required_filtered.fasta"
    table = OUT / "all_run_counts.tsv"
    if not merged.exists():
        run(["vsearch", "--fastq_mergepairs", str(forward), "--reverse", str(reverse),
             "--fastqout", str(merged), "--fastq_minovlen", "16", "--fastq_maxdiffs", "300",
             "--threads", "16"])
    if not filtered.exists():
        run(["vsearch", "--fastq_filter", str(merged), "--fastq_maxee", "0.1",
             "--fastq_minlen", "100", "--fastaout", str(filtered)])
    if not table.exists():
        map_with_minimap2(filtered, table, rows)


if __name__ == "__main__":
    main()
