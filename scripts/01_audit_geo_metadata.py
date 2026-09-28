#!/usr/bin/env python3
"""Audit GEO series-matrix metadata for the two manuscript datasets.

Run: python3 scripts/01_audit_geo_metadata.py --cache-dir <download-folder>
The cache directory is optional. This script uses only the Python standard library.
It never infers animal identity, tissue side, or processing batch from replicate labels.
"""

import argparse
import csv
import gzip
import hashlib
import re
import urllib.request
from collections import Counter
from pathlib import Path


ACCESSIONS = ("GSE160543", "GSE246156")
FIELDS = ("accession", "gsm", "title", "condition", "tissue", "drg_level",
          "day", "replicate_label", "animal_id", "side", "sra_accession")


def source_url(accession):
    prefix = accession[:-3] + "nnn"
    return (f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/"
            f"{accession}/matrix/{accession}_series_matrix.txt.gz")


def read_bytes(accession, cache_dir):
    path = cache_dir / f"{accession}_series_matrix.txt.gz"
    if path.exists():
        data = path.read_bytes()
    else:
        with urllib.request.urlopen(source_url(accession), timeout=60) as response:
            data = response.read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    return data


def parse_matrix(data):
    rows = {}
    for line in gzip.decompress(data).decode("utf-8").splitlines():
        if line == "!series_matrix_table_begin":
            break
        if line.startswith("!Sample_"):
            parts = next(csv.reader([line], delimiter="\t"))
            rows[parts[0]] = parts[1:]
    count = len(rows["!Sample_geo_accession"])
    if not all(len(values) == count for values in rows.values()):
        raise ValueError("Inconsistent GEO sample metadata columns")
    return rows


def sample_record(accession, rows, index):
    title = rows["!Sample_title"][index]
    gsm = rows["!Sample_geo_accession"][index]
    relation = rows["!Sample_relation"][index]
    sra = re.search(r"SRX\d+", relation)
    if accession == "GSE160543":
        match = re.fullmatch(r"(Vehicle|Paclitaxel|Oxaliplatin) rep(\d+)",
                             title, flags=re.IGNORECASE)
        if not match:
            raise ValueError(f"Unrecognized {accession} title: {title}")
        condition, replicate = match.groups()
        condition = condition.capitalize()
        tissue, level, day = "DRG", "", ""
    else:
        match = re.fullmatch(
            r"(L[456] DRGs|Sciatic nerve),"
            r"(L5 Nerve entrapment|control) for ([37]) days,rep(\d+)", title)
        if not match:
            raise ValueError(f"Unrecognized {accession} title: {title}")
        label, intervention, day, replicate = match.groups()
        condition = "Compression" if intervention == "L5 Nerve entrapment" else "Sham"
        tissue = "Sciatic nerve" if label == "Sciatic nerve" else "DRG"
        level = label.split()[0] if tissue == "DRG" else ""
    return dict(accession=accession, gsm=gsm, title=title, condition=condition,
                tissue=tissue, drg_level=level, day=day,
                replicate_label=replicate, animal_id="", side="",
                sra_accession=sra.group() if sra else "")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=Path("data/source-cache"))
    parser.add_argument("--out-dir", type=Path, default=Path("data/metadata"))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    provenance = []
    records = []
    for accession in ACCESSIONS:
        data = read_bytes(accession, args.cache_dir)
        rows = parse_matrix(data)
        records.extend(sample_record(accession, rows, i)
                       for i in range(len(rows["!Sample_geo_accession"])))
        provenance.append((accession, source_url(accession),
                           hashlib.sha256(data).hexdigest(), len(data)))
    assert len({r["gsm"] for r in records}) == len(records)
    assert Counter(r["accession"] for r in records) == {
        "GSE160543": 12, "GSE246156": 48}
    with (args.out_dir / "geo_sample_manifest.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(records)
    with (args.out_dir / "source_provenance.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("accession", "source_url", "sha256_gz", "size_bytes"))
        writer.writerows(provenance)
    counts = Counter((r["accession"], r["condition"], r["tissue"],
                      r["drg_level"], r["day"]) for r in records)
    with (args.out_dir / "sample_groups.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(("accession", "condition", "tissue", "drg_level", "day", "n_samples"))
        for key, n in sorted(counts.items()):
            writer.writerow((*key, n))
    print(f"Wrote {len(records)} sample rows and {len(counts)} design cells")


if __name__ == "__main__":
    main()
