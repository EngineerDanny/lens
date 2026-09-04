#!/usr/bin/env python3
"""Download Carlström et al. 2019 reads and accompanying analysis files."""

from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import time
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "external_data" / "carlstrom_phyllosphere_2019"
READS = OUT / "raw_reads"
ENA_REPORT = OUT / "ena_read_runs.tsv"
REPOSITORY = OUT / "source"
ENA_URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB32997"
    "&result=read_run&fields=run_accession,run_alias,experiment_alias,sample_alias,"
    "sample_title,description,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&download=true"
)


def retrieve(url: str, destination: Path) -> None:
    temporary = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as response, temporary.open("wb") as handle:
                shutil.copyfileobj(response, handle, length=1024 * 1024)
            temporary.replace(destination)
            return
        except Exception:
            temporary.unlink(missing_ok=True)
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_one(item):
    url, destination, expected_md5, expected_size = item
    if destination.exists() and destination.stat().st_size == expected_size:
        if md5sum(destination) == expected_md5:
            return destination.name, "present"
    retrieve(url, destination)
    if destination.stat().st_size != expected_size:
        raise RuntimeError(f"Size mismatch for {destination.name}")
    if md5sum(destination) != expected_md5:
        raise RuntimeError(f"MD5 mismatch for {destination.name}")
    return destination.name, "downloaded"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    READS.mkdir(parents=True, exist_ok=True)
    if not REPOSITORY.exists():
        subprocess.run(
            ["git", "clone", "--depth", "1", "https://github.com/cmfield/carlstrom2019.git", str(REPOSITORY)],
            check=True,
        )
    retrieve(ENA_URL, ENA_REPORT)
    with ENA_REPORT.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    jobs = []
    for row in rows:
        urls = ["https://" + value for value in row["fastq_ftp"].split(";")]
        checksums = row["fastq_md5"].split(";")
        sizes = [int(value) for value in row["fastq_bytes"].split(";")]
        for index, (url, checksum, size) in enumerate(zip(urls, checksums, sizes), start=1):
            destination = READS / f"{row['run_accession']}_{index}.fastq.gz"
            jobs.append((url, destination, checksum, size))
    completed = 0
    with ThreadPoolExecutor(max_workers=12) as executor:
        futures = [executor.submit(download_one, job) for job in jobs]
        for future in as_completed(futures):
            future.result()
            completed += 1
            if completed % 100 == 0 or completed == len(jobs):
                print(f"Verified {completed}/{len(jobs)} files", flush=True)


if __name__ == "__main__":
    main()
