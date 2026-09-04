#!/usr/bin/env python3
"""Map Carlström et al. reads to the 62 strain reference sequences."""

from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "external_data" / "carlstrom_phyllosphere_2019"
READS = BASE / "raw_reads"
REPORT = BASE / "ena_read_runs.tsv"
REFERENCE = BASE / "source" / "64_strains_sequences.fasta"
OUT = BASE / "processed"
RUN_COUNTS = OUT / "run_counts"


def process(row):
    accession = row["run_accession"]
    destination = RUN_COUNTS / f"{accession}.tsv"
    if destination.exists() and destination.stat().st_size > 0:
        return accession, "present"
    with tempfile.TemporaryDirectory(prefix=f"{accession}.", dir=OUT / "tmp") as temporary:
        temporary = Path(temporary)
        merged = temporary / "merged.fastq"
        filtered = temporary / "filtered.fasta"
        table = temporary / "counts.tsv"
        subprocess.run([
            "vsearch", "--fastq_mergepairs", str(READS / f"{accession}_1.fastq.gz"),
            "--reverse", str(READS / f"{accession}_2.fastq.gz"),
            "--fastqout", str(merged), "--fastq_minovlen", "16",
            "--fastq_maxdiffs", "300", "--threads", "1", "--quiet",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            "vsearch", "--fastq_filter", str(merged), "--fastq_maxee", "0.1",
            "--fastq_minlen", "100", "--fastaout", str(filtered), "--quiet",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([
            "vsearch", "--usearch_global", str(filtered), "--db", str(REFERENCE),
            "--id", "0.97", "--strand", "both", "--top_hits_only",
            "--otutabout", str(table), "--threads", "1", "--quiet",
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        table.replace(destination)
    return accession, "processed"


def main():
    RUN_COUNTS.mkdir(parents=True, exist_ok=True)
    (OUT / "tmp").mkdir(parents=True, exist_ok=True)
    with REPORT.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required_titles = set()
    for experiment in (53, 61, 62):
        with (BASE / "source" / f"metadata{experiment}.csv").open() as handle:
            metadata = list(csv.DictReader(handle, delimiter=";"))
        for record in metadata:
            if experiment == 53:
                keep = record["Treatment"] == "ALL" and record["Spray"] in {"U", "Mg"} and record["Time"] != "t0"
            else:
                keep = record["Treatment"] not in {"Ax", "Axenic", "CORE"} and record["Time"] != "t0"
            if keep:
                required_titles.add(record["Name"])
    rows = [row for row in rows if row["sample_title"] in required_titles]
    rows = [
        row for row in rows
        if (READS / f"{row['run_accession']}_1.fastq.gz").exists()
        and (READS / f"{row['run_accession']}_2.fastq.gz").exists()
    ]
    completed = 0
    with ThreadPoolExecutor(max_workers=16) as executor:
        futures = [executor.submit(process, row) for row in rows]
        for future in as_completed(futures):
            future.result()
            completed += 1
            if completed % 50 == 0 or completed == len(rows):
                print(f"Processed {completed}/{len(rows)} runs", flush=True)


if __name__ == "__main__":
    main()
