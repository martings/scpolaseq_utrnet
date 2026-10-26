"""Lossless input keys and schema validation."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

CONTRACT_VERSION = "utrnet-input/v0.1"
REQUIRED = {
    "cells": {"sample_id", "library_id", "barcode_corrected", "cell_id"},
    "pas_sites": {"site_id", "gene_id", "chrom", "start", "end", "strand"},
    "pas_counts": {"library_id", "group_level", "group_id", "site_id", "cell_barcode", "umi_count"},
}


def require_columns(frame: pd.DataFrame, kind: str) -> None:
    missing = REQUIRED[kind] - set(frame.columns)
    if missing:
        raise ValueError(f"{kind}: missing columns {sorted(missing)}")


def read_tsv(path: Path, *, usecols: list[str] | None = None) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=usecols)


def validate_cells(cells: pd.DataFrame) -> None:
    require_columns(cells, "cells")
    key = ["sample_id", "library_id", "barcode_corrected"]
    if cells.duplicated(key).any() or cells.cell_id.duplicated().any():
        raise ValueError("duplicate exact cell key or cell_id")
    expected = cells.sample_id + ":" + cells.library_id + ":" + cells.barcode_corrected
    if not expected.equals(cells.cell_id):
        raise ValueError("cell_id is not reversible from raw sample/library/barcode")


def validate_counts(counts: pd.DataFrame, cells: pd.DataFrame, sites: pd.DataFrame) -> None:
    require_columns(counts, "pas_counts")
    if not bool(counts.site_id.isin(sites.site_id).all()):
        raise ValueError("Sierra site_id absent from PAS reference")
    key = cells[["library_id", "barcode_corrected"]]
    if key.duplicated().any():
        raise ValueError("library/barcode collision across samples")
    check = counts[["library_id", "cell_barcode"]].drop_duplicates().merge(
        key, left_on=["library_id", "cell_barcode"],
        right_on=["library_id", "barcode_corrected"], how="left", indicator=True,
    )
    if (check._merge != "both").any():
        raise ValueError("Sierra cell barcode absent from exact library/barcode registry")
    numeric = pd.to_numeric(counts.umi_count, errors="raise")
    if numeric.isna().any() or (numeric < 0).any() or (numeric % 1 != 0).any():
        raise ValueError("invalid UMI count")
