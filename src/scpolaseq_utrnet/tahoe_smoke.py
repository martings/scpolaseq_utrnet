"""End-to-end structural smoke test for the checked-in synthetic fixture."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from .tahoe import sparse_counts


def run_fixture(fixture_path: Path) -> dict[str, Any]:
    fixture = json.loads(Path(fixture_path).read_text(encoding="utf-8"))
    gene_metadata = pd.DataFrame(fixture["gene_metadata"])
    matrix, cell_ids, gene_ids = sparse_counts(
        fixture["expression_records"],
        gene_metadata,
    )
    expected_count = sum(
        int(value) for row in fixture["expression_records"] for value in row["expressions"][1:]
    )
    observed_count = int(matrix.sum())
    if expected_count != observed_count:
        raise ValueError("smoke fixture count conservation failed")
    return {
        "status": "pass",
        "purpose": fixture["purpose"],
        "source_dataset": fixture["dataset_id"],
        "source_revision": fixture["source_revision"],
        "cells": len(cell_ids),
        "genes": len(gene_ids),
        "nonzero_counts": int(matrix.nnz),
        "raw_count_sum_after_marker_exclusion": observed_count,
        "missing_modalities": ["PAS usage", "timepoint"],
        "biological_claim": False,
    }


def main() -> int:
    fixture = Path(__file__).resolve().parents[2] / "examples" / "tahoe_smoke_fixture.json"
    print(json.dumps(run_fixture(fixture), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
