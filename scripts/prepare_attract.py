"""Verify the ATtRACT 2016 archive and regenerate the two-row MVP motif subset."""
from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
import zipfile


ARCHIVE_SHA256 = "adf73651e8e4d75c0fea86a5fb10926c93b7b04628d065789c1196b9ee0a3a76"
SOURCE_URL = "https://attract.cnic.es/download"
TARGETS = {
    ("ELAVL1", "AUUUA", "9882309"): ("ATtRACT_1103_ELAVL1_AUUUA", "HeLa nuclear extracts and recombinant protein"),
    ("ZFP36", "AUUUA", "18256032"): ("ATtRACT_38_ZFP36_AUUUA", "other or unspecified context"),
}


def prepare(archive: Path, output: Path) -> None:
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != ARCHIVE_SHA256:
        raise ValueError(f"ATtRACT archive SHA-256 mismatch: {digest}")
    selected = {}
    with zipfile.ZipFile(archive) as bundle:
        with bundle.open("ATtRACT_db.txt") as stream:
            rows = csv.DictReader((line.decode() for line in stream), delimiter="\t")
            for row in rows:
                key = (row["Gene_name"], row["Motif"], row["Pubmed"])
                if key in TARGETS and row["Organism"] == "Homo_sapiens":
                    selected[key] = row
    if set(selected) != set(TARGETS):
        raise ValueError(f"missing curated ATtRACT records: {set(TARGETS) - set(selected)}")
    columns = ["motif_id", "regulator", "regulator_class", "pattern", "evidence_tier", "source",
               "version", "url", "pubmed_id", "assay", "cell_context", "matrix_id"]
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in TARGETS:
            row = selected[key]
            motif_id, context = TARGETS[key]
            writer.writerow({"motif_id": motif_id, "regulator": row["Gene_name"],
                             "regulator_class": "ARE-BP", "pattern": row["Motif"],
                             "evidence_tier": "experimental_other_context", "source": "ATtRACT",
                             "version": "2016-09-27", "url": SOURCE_URL, "pubmed_id": row["Pubmed"],
                             "assay": row["Experiment_description"], "cell_context": context,
                             "matrix_id": row["Matrix_id"]})


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    prepare(args.archive, args.output)
