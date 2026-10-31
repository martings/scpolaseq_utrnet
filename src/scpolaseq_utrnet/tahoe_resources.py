"""Measured sparse payload reporting and transparent GNN memory estimates."""
from __future__ import annotations

from typing import Any

from scipy import sparse


def sparse_resource_report(
    matrix: sparse.spmatrix,
    *,
    elapsed_seconds: float | None = None,
    observed_peak_ram_bytes: int | None = None,
    observed_peak_vram_bytes: int | None = None,
    cache_bytes: int | None = None,
) -> dict[str, Any]:
    """Report measured CSR storage separately from optional run observations."""
    csr = matrix.tocsr()
    payload_bytes = int(csr.data.nbytes + csr.indices.nbytes + csr.indptr.nbytes)
    result: dict[str, Any] = {
        "shape": [int(csr.shape[0]), int(csr.shape[1])],
        "nnz": int(csr.nnz),
        "csr_payload_bytes": payload_bytes,
        "csr_payload_is_ram_peak": False,
        "observed_peak_ram_bytes": observed_peak_ram_bytes,
        "observed_peak_vram_bytes": observed_peak_vram_bytes,
        "cache_bytes": cache_bytes,
        "elapsed_seconds": elapsed_seconds,
        "cells_per_second": (
            float(csr.shape[0] / elapsed_seconds)
            if elapsed_seconds is not None and elapsed_seconds > 0 else None
        ),
        "genes_per_second": (
            float(csr.nnz / elapsed_seconds)
            if elapsed_seconds is not None and elapsed_seconds > 0 else None
        ),
    }
    return result


def estimate_gnn_activation_bytes(
    *,
    batch_size: int,
    n_genes: int,
    latent_dim: int = 64,
    layers: int = 2,
    bytes_per_value: int = 4,
) -> dict[str, int]:
    """Upper-bound node activations for a dense per-sample message-passing batch."""
    if min(batch_size, n_genes, latent_dim, layers, bytes_per_value) < 1:
        raise ValueError("estimate dimensions must be positive")
    node_tensor_count = layers + 2
    activation_bytes = batch_size * n_genes * latent_dim * bytes_per_value * node_tensor_count
    return {
        "batch_size": batch_size,
        "n_genes": n_genes,
        "latent_dim": latent_dim,
        "layers": layers,
        "node_tensor_count_estimate": node_tensor_count,
        "activation_bytes_estimate": activation_bytes,
    }
