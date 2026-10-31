"""Loss-aware Tahoe-100M metadata and sparse-count adapter.

The module never downloads data at import time. Remote scans must be explicitly
started by a caller and should use a pinned dataset revision and a byte budget.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
import ast
import json
import re
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse

DATASET_ID = "tahoebio/Tahoe-100M"
DATASET_REVISION = "349e2a55f5538a3e73b32a06c6bdd8e9e87011e8"
CONFIGS = (
    "expression_data",
    "sample_metadata",
    "gene_metadata",
    "drug_metadata",
    "cell_line_metadata",
    "obs_metadata",
    "pseudobulk_differential_expression",
)
EXPRESSION_COLUMNS = {
    "genes", "expressions", "drug", "sample", "BARCODE_SUB_LIB_ID",
    "cell_line_id", "moa-fine", "canonical_smiles", "pubchem_cid", "plate",
}
SAMPLE_COLUMNS = {"sample", "plate", "drug", "drugname_drugconc"}
GENE_COLUMNS = {"token_id", "gene_symbol", "ensembl_id"}
CONTROL_DRUG = "DMSO_TF"


def validate_columns(frame: pd.DataFrame, required: set[str], table: str) -> None:
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{table}: missing required columns {sorted(missing)}")


def make_gene_map(gene_metadata: pd.DataFrame) -> dict[int, str]:
    """Map Tahoe tokens to stable Ensembl gene IDs without symbol substitution."""
    validate_columns(gene_metadata, GENE_COLUMNS, "gene_metadata")
    mapping: dict[int, str] = {}
    for token, ensembl in zip(gene_metadata.token_id, gene_metadata.ensembl_id):
        token_id = int(token)
        gene_id = str(ensembl)
        if not gene_id or gene_id.lower() == "nan":
            continue
        if token_id in mapping and mapping[token_id] != gene_id:
            raise ValueError(f"token_id {token_id} maps to multiple Ensembl IDs")
        mapping[token_id] = gene_id
    if not mapping:
        raise ValueError("gene_metadata produced an empty token map")
    return mapping


def parse_concentrations(value: Any) -> list[dict[str, Any]]:
    """Parse the documented drugname_drugconc list while retaining units."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    if isinstance(value, list):
        parsed = value
    else:
        try:
            parsed = ast.literal_eval(str(value))
        except (ValueError, SyntaxError):
            return [{"raw": str(value)}]
    if not isinstance(parsed, (list, tuple)):
        return [{"raw": str(value)}]
    output = []
    for entry in parsed:
        if isinstance(entry, (list, tuple)) and len(entry) >= 3:
            output.append({"drug": str(entry[0]), "concentration": entry[1], "unit": str(entry[2])})
        else:
            output.append({"raw": str(entry)})
    return output


def expression_row_to_pairs(
    row: Mapping[str, Any],
    token_to_ensembl: Mapping[int, str],
    *,
    marker_token_id: int | None = None,
) -> list[tuple[str, int]]:
    """Drop exactly the leading marker pair, then validate and map raw counts."""
    tokens = list(row["genes"])
    values = list(row["expressions"])
    if len(tokens) != len(values):
        raise ValueError("genes and expressions lengths differ")
    if not tokens:
        raise ValueError("expression row is empty")
    leading = int(tokens[0])
    if marker_token_id is not None and leading != marker_token_id:
        raise ValueError(f"unexpected marker token {leading}; expected {marker_token_id}")
    pairs: list[tuple[str, int]] = []
    seen: set[str] = set()
    for token, value in zip(tokens[1:], values[1:]):
        token_id = int(token)
        if token_id not in token_to_ensembl:
            raise ValueError(f"unmapped gene token {token_id}")
        count = float(value)
        if not np.isfinite(count) or count < 0 or not count.is_integer():
            raise ValueError(f"invalid raw count for token {token_id}: {value}")
        gene_id = token_to_ensembl[token_id]
        if gene_id in seen:
            raise ValueError(f"duplicate gene ID {gene_id} in one cell")
        seen.add(gene_id)
        pairs.append((gene_id, int(count)))
    return pairs


def sparse_counts(
    rows: Iterable[Mapping[str, Any]],
    gene_metadata: pd.DataFrame,
    *,
    marker_token_id: int | None = None,
) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    """Build a cell-by-gene CSR matrix; row order and raw cell IDs are preserved."""
    token_map = make_gene_map(gene_metadata)
    gene_ids = sorted(set(token_map.values()))
    gene_col = {gene_id: i for i, gene_id in enumerate(gene_ids)}
    cell_ids: list[str] = []
    row_ix: list[int] = []
    col_ix: list[int] = []
    values: list[int] = []
    seen_cells: set[str] = set()
    for row in rows:
        cell_id = str(row["BARCODE_SUB_LIB_ID"])
        if cell_id in seen_cells:
            raise ValueError(f"duplicate BARCODE_SUB_LIB_ID {cell_id}")
        seen_cells.add(cell_id)
        cell_ix = len(cell_ids)
        cell_ids.append(cell_id)
        for gene_id, count in expression_row_to_pairs(row, token_map, marker_token_id=marker_token_id):
            row_ix.append(cell_ix)
            col_ix.append(gene_col[gene_id])
            values.append(count)
    matrix = sparse.coo_matrix(
        (np.asarray(values, dtype=np.int64), (row_ix, col_ix)),
        shape=(len(cell_ids), len(gene_ids)),
        dtype=np.int64,
    ).tocsr()
    return matrix, cell_ids, gene_ids


def select_metadata_first_pilot(
    sample_index: pd.DataFrame,
    *,
    cell_lines: int = 2,
    conditions: int = 4,
    min_treatment_samples: int = 2,
    min_control_samples: int = 2,
    min_cells: int = 10,
) -> dict[str, Any]:
    """Select replicated conditions only when plate-and-cell-line controls exist.

    sample_index is a small, metadata-only table with one row per sample, built
    from dataset metadata. Sample count is the replicate unit; cell count is a
    quality/coverage field and is never treated as biological replication.
    """
    if cell_lines < 1 or not 3 <= conditions <= 5:
        raise ValueError("pilot requires at least one cell line and 3–5 conditions")
    required = {"sample", "plate", "cell_line_id", "drug", "drugname_drugconc", "n_cells"}
    validate_columns(sample_index, required, "sample_index")
    if sample_index["sample"].astype(str).duplicated().any():
        raise ValueError("sample_index must have one unique row per sample")
    frame = sample_index.copy()
    for column in ("sample", "plate", "cell_line_id", "drug", "drugname_drugconc"):
        frame[column] = frame[column].astype(str)
    frame["n_cells"] = pd.to_numeric(frame["n_cells"], errors="raise")
    controls = frame[frame.drug == CONTROL_DRUG]
    treatments = frame[frame.drug != CONTROL_DRUG]
    candidates: list[dict[str, Any]] = []
    group_cols = ["cell_line_id", "drug", "drugname_drugconc"]
    for key, group in treatments.groupby(group_cols, sort=True):
        cell_line, drug, dose_label = key
        group = group[group.n_cells >= min_cells]
        matched = controls[controls.cell_line_id == cell_line]
        treatment_plates = sorted(set(group.plate))
        matched = matched[matched.plate.isin(treatment_plates)]
        if group["sample"].nunique() < min_treatment_samples or matched["sample"].nunique() < min_control_samples:
            continue
        treatment_samples = sorted(group["sample"].unique())
        control_samples = sorted(matched["sample"].unique())
        if len(treatment_samples) < min_treatment_samples:
            continue
        candidates.append({
            "cell_line_id": cell_line,
            "drug": drug,
            "condition": dose_label,
            "concentrations": parse_concentrations(dose_label),
            "treatment_samples": treatment_samples,
            "control_samples": control_samples,
            "matched_plates": treatment_plates,
            "treatment_sample_n": len(treatment_samples),
            "control_sample_n": len(control_samples),
            "treatment_cell_n": int(group.n_cells.sum()),
            "control_cell_n": int(matched.n_cells.sum()),
            "selection_basis": "unique sample IDs are replicate units; cell counts are coverage only",
        })
    candidates.sort(key=lambda x: (
        -min(x["treatment_sample_n"], x["control_sample_n"]),
        -x["treatment_cell_n"],
        x["cell_line_id"], x["drug"], x["condition"],
    ))
    selected: list[dict[str, Any]] = []
    used_lines: set[str] = set()
    for item in candidates:
        if item["cell_line_id"] not in used_lines and len(used_lines) < cell_lines:
            selected.append(item)
            used_lines.add(item["cell_line_id"])
    if len(used_lines) < cell_lines:
        raise ValueError(f"only {len(used_lines)} cell lines meet matched-control/replicate criteria")
    for item in candidates:
        if len(selected) >= conditions:
            break
        if item not in selected and item["cell_line_id"] in used_lines:
            selected.append(item)
    if len(selected) < 3:
        raise ValueError(f"only {len(selected)} replicated, plate-matched conditions qualify")
    return {
        "contract": "cell-regulatory-state/v0.1",
        "source_dataset": DATASET_ID,
        "source_revision": DATASET_REVISION,
        "selection_method": "metadata-first, deterministic quality/replicate ranking",
        "sample_is_replicate_unit": True,
        "cell_count_is_biological_replication": False,
        "conditions": selected,
    }


def validate_group_partition(
    assignments: pd.DataFrame,
    *,
    group_column: str,
    partition_column: str = "partition",
) -> None:
    """Reject any sample/replicate group that appears in more than one split."""
    validate_columns(assignments, {group_column, partition_column}, "partition assignments")
    leaks = assignments.groupby(group_column)[partition_column].nunique()
    if (leaks > 1).any():
        raise ValueError(f"partition leakage for groups {list(leaks[leaks > 1].index)}")
