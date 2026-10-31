"""Small, explicit-budget sparse pilot cache writer."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Mapping

from scipy import sparse

from .tahoe import DATASET_ID, DATASET_REVISION


def write_sparse_cache(
    matrix: sparse.csr_matrix,
    cell_ids: list[str],
    gene_ids: list[str],
    output: Path,
    *,
    max_bytes: int,
    source_revision: str = DATASET_REVISION,
    metadata: Mapping[str, Any] | None = None,
) -> tuple[Path, Path]:
    """Write a sparse matrix plus identity sidecar under an explicit byte cap."""
    if max_bytes <= 0:
        raise ValueError("max_bytes must be explicitly positive")
    if matrix.shape != (len(cell_ids), len(gene_ids)):
        raise ValueError("matrix dimensions do not match the supplied identifiers")
    if len(set(cell_ids)) != len(cell_ids) or len(set(gene_ids)) != len(gene_ids):
        raise ValueError("cell and gene identifiers must be unique")
    target = Path(output)
    target.parent.mkdir(parents=True, exist_ok=True)
    matrix_tmp = target.with_name(target.name + ".tmp.npz")
    sidecar = target.with_name(target.name + ".json")
    sidecar_tmp = sidecar.with_name(sidecar.name + ".tmp")
    sparse.save_npz(matrix_tmp, matrix.tocsr())
    record = {
        "contract": "cell-regulatory-state/v0.1",
        "source_dataset": DATASET_ID,
        "source_revision": source_revision,
        "matrix_path": target.name,
        "shape": list(matrix.shape),
        "cell_ids": cell_ids,
        "gene_ids": gene_ids,
        "raw_count_sum": int(matrix.sum()),
        "nonzero_count": int(matrix.nnz),
        "missing_modalities": ["PAS usage", "timepoint"],
        "metadata": dict(metadata or {}),
    }
    sidecar_tmp.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    size = matrix_tmp.stat().st_size + sidecar_tmp.stat().st_size
    if size > max_bytes:
        matrix_tmp.unlink(missing_ok=True)
        sidecar_tmp.unlink(missing_ok=True)
        raise ValueError(f"pilot cache would use {size} bytes, exceeding max_bytes={max_bytes}")
    matrix_tmp.replace(target)
    sidecar_tmp.replace(sidecar)
    return target, sidecar
