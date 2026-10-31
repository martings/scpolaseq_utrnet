"""Pinned, metadata-only Hugging Face inspection helpers for Tahoe-100M."""
from __future__ import annotations

from collections import Counter
from typing import Any, Iterable, Mapping

import pandas as pd

from .tahoe import CONFIGS, DATASET_ID, DATASET_REVISION, SAMPLE_COLUMNS, validate_columns


def inspect_dataset_revision(revision: str = DATASET_REVISION) -> dict[str, Any]:
    """Read repo and builder metadata only; never fetch a Parquet data shard."""
    try:
        from datasets import get_dataset_config_names, load_dataset_builder
        from huggingface_hub import HfApi
    except ImportError as exc:
        raise RuntimeError("install the optional tahoe extra to inspect Hugging Face metadata") from exc

    api = HfApi()
    info = api.dataset_info(DATASET_ID, revision=revision, files_metadata=True)
    resolved_revision = getattr(info, "sha", None) or revision
    siblings = getattr(info, "siblings", None) or []
    files = []
    for sibling in siblings:
        path = getattr(sibling, "rfilename", None)
        size = getattr(sibling, "size", None)
        if path:
            files.append({"path": path, "size_bytes": size})
    config_names = get_dataset_config_names(DATASET_ID, revision=resolved_revision)
    schemas: dict[str, Any] = {}
    for config in config_names:
        builder = load_dataset_builder(DATASET_ID, config, revision=resolved_revision)
        features = builder.info.features
        schemas[config] = features.to_dict() if features is not None else None
    return {
        "dataset_id": DATASET_ID,
        "requested_revision": revision,
        "resolved_revision": resolved_revision,
        "configs": list(config_names),
        "expected_configs": list(CONFIGS),
        "missing_expected_configs": sorted(set(CONFIGS) - set(config_names)),
        "schemas": schemas,
        "files": files,
        "data_shards_downloaded": 0,
    }


def sample_index_from_metadata_rows(
    expression_metadata_rows: Iterable[Mapping[str, Any]],
    sample_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate projected cell metadata per source sample and cell line.

    Tahoe samples contain several co-cultured cell lines, so sample alone is
    not a row key. Replicate counts remain unique sample IDs within each line.
    Both expression_data.cell_line_id and obs_metadata.cell_line are accepted.
    """
    validate_columns(sample_metadata, SAMPLE_COLUMNS, "sample_metadata")
    samples = sample_metadata.copy()
    samples["sample"] = samples["sample"].astype(str)
    if samples["sample"].duplicated().any():
        raise ValueError("sample_metadata has duplicate sample keys")
    metadata = samples.set_index("sample", drop=False)
    cell_counts: Counter[tuple[str, str]] = Counter()
    observed: dict[tuple[str, str], dict[str, set[str]]] = {}
    seen_cell_ids: set[str] = set()
    qc: dict[tuple[str, str], dict[str, float]] = {}
    for row in expression_metadata_rows:
        sample = str(row["sample"])
        cell_id = str(row["BARCODE_SUB_LIB_ID"])
        cell_line = str(row.get("cell_line_id", row.get("cell_line", "")))
        if not cell_line:
            raise ValueError("metadata row lacks cell_line_id/cell_line")
        if cell_id in seen_cell_ids:
            raise ValueError(f"duplicate BARCODE_SUB_LIB_ID in metadata: {cell_id}")
        seen_cell_ids.add(cell_id)
        if sample not in metadata.index:
            raise ValueError(f"metadata sample absent from sample_metadata: {sample}")
        key = (sample, cell_line)
        cell_counts[key] += 1
        state = observed.setdefault(key, {"plate": set(), "drug": set(), "drugname_drugconc": set()})
        for column in ("plate", "drug"):
            state[column].add(str(row[column]))
        if "drugname_drugconc" in row:
            state["drugname_drugconc"].add(str(row["drugname_drugconc"]))
        q = qc.setdefault(key, {"gene_sum": 0.0, "gene_n": 0.0, "mito_sum": 0.0, "mito_n": 0.0})
        gene_count = row.get("gene_count", row.get("mean_gene_count"))
        mito = row.get("pcnt_mito", row.get("mean_pcnt_mito"))
        if gene_count is not None:
            q["gene_sum"] += float(gene_count)
            q["gene_n"] += 1
        if mito is not None:
            q["mito_sum"] += float(mito)
            q["mito_n"] += 1
    records = []
    for (sample, cell_line), state in observed.items():
        if any(len(values) != 1 for values in state.values()):
            raise ValueError(f"inconsistent metadata within sample {sample}, cell line {cell_line}")
        source = metadata.loc[sample]
        plate = next(iter(state["plate"]))
        drug = next(iter(state["drug"]))
        if plate != str(source["plate"]) or drug != str(source["drug"]):
            raise ValueError(f"sample_metadata conflict for sample {sample}")
        if state["drugname_drugconc"] and next(iter(state["drugname_drugconc"])) != str(source["drugname_drugconc"]):
            raise ValueError(f"dose metadata conflict for sample {sample}")
        q = qc[(sample, cell_line)]
        records.append({
            "sample": sample,
            "plate": plate,
            "cell_line_id": cell_line,
            "drug": drug,
            "drugname_drugconc": str(source["drugname_drugconc"]),
            "n_cells": int(cell_counts[(sample, cell_line)]),
            "mean_gene_count": q["gene_sum"] / q["gene_n"] if q["gene_n"] else source.get("mean_gene_count"),
            "mean_tscp_count": source.get("mean_tscp_count"),
            "mean_pcnt_mito": q["mito_sum"] / q["mito_n"] if q["mito_n"] else source.get("mean_pcnt_mito"),
        })
    return pd.DataFrame.from_records(records)
def inspect_file_budget(audit: Mapping[str, Any], *, max_metadata_bytes: int) -> dict[str, Any]:
    """Summarize the known metadata-file cache ceiling before a caller downloads it."""
    if max_metadata_bytes <= 0:
        raise ValueError("max_metadata_bytes must be explicitly positive")
    metadata_files = [
        item for item in audit.get("files", [])
        if str(item.get("path", "")).startswith("metadata/") and not str(item.get("path", "")).startswith("metadata/pseudobulk_differential_expression/")
    ]
    known_bytes = sum(int(item["size_bytes"]) for item in metadata_files if item.get("size_bytes") is not None)
    return {
        "max_metadata_bytes": max_metadata_bytes,
        "known_metadata_bytes": known_bytes,
        "within_budget": known_bytes <= max_metadata_bytes,
        "metadata_file_count": len(metadata_files),
        "expression_bytes_downloaded": 0,
        "pseudobulk_de_bytes_downloaded": 0,
    }
