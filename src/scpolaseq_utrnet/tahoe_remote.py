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
    metadata_rows: Iterable[Mapping[str, Any]],
    sample_metadata: pd.DataFrame,
) -> pd.DataFrame:
    """Aggregate projected cell metadata per source sample and cell line.

    Tahoe samples contain several co-cultured cell lines, so sample alone is
    not a row key. The expression table contains full-filter cells only;
    obs_metadata also contains minimal-filter cells. n_cells therefore counts
    only full rows while the additional columns preserve both source totals.
    """
    validate_columns(sample_metadata, SAMPLE_COLUMNS, "sample_metadata")
    samples = sample_metadata.copy()
    samples["sample"] = samples["sample"].astype(str)
    if samples["sample"].duplicated().any():
        raise ValueError("sample_metadata has duplicate sample keys")
    metadata = samples.set_index("sample", drop=False)
    totals: dict[tuple[str, str], dict[str, Any]] = {}
    for row in metadata_rows:
        sample = str(row["sample"])
        cell_line = str(row.get("cell_line_id", row.get("cell_line", "")))
        if not cell_line:
            raise ValueError("metadata row lacks cell_line_id/cell_line")
        if sample not in metadata.index:
            raise ValueError(f"metadata sample absent from sample_metadata: {sample}")
        key = (sample, cell_line)
        state = totals.setdefault(key, {
            "plate": set(), "drug": set(), "drugname_drugconc": set(),
            "n_cells": 0, "n_cells_minimal": 0, "n_cells_obs_total": 0,
            "gene_sum": 0.0, "gene_n": 0, "mito_sum": 0.0, "mito_n": 0,
        })
        state["plate"].add(str(row["plate"]))
        state["drug"].add(str(row["drug"]))
        if "drugname_drugconc" in row:
            state["drugname_drugconc"].add(str(row["drugname_drugconc"]))
        pass_filter = str(row.get("pass_filter", "full")).lower()
        state["n_cells_obs_total"] += 1
        if pass_filter == "minimal":
            state["n_cells_minimal"] += 1
        elif pass_filter == "full":
            state["n_cells"] += 1
            gene_count = row.get("gene_count", row.get("mean_gene_count"))
            mito = row.get("pcnt_mito", row.get("mean_pcnt_mito"))
            if gene_count is not None:
                state["gene_sum"] += float(gene_count)
                state["gene_n"] += 1
            if mito is not None:
                state["mito_sum"] += float(mito)
                state["mito_n"] += 1
        else:
            raise ValueError(f"unknown obs_metadata.pass_filter value: {pass_filter}")
    records = []
    for (sample, cell_line), state in totals.items():
        if any(len(state[column]) != 1 for column in ("plate", "drug")):
            raise ValueError(f"inconsistent metadata within sample {sample}, cell line {cell_line}")
        if state["drugname_drugconc"] and len(state["drugname_drugconc"]) != 1:
            raise ValueError(f"inconsistent dose metadata within sample {sample}, cell line {cell_line}")
        source = metadata.loc[sample]
        plate = next(iter(state["plate"]))
        drug = next(iter(state["drug"]))
        if plate != str(source["plate"]) or drug != str(source["drug"]):
            raise ValueError(f"sample_metadata conflict for sample {sample}")
        dose = str(source["drugname_drugconc"])
        if state["drugname_drugconc"] and next(iter(state["drugname_drugconc"])) != dose:
            raise ValueError(f"dose metadata conflict for sample {sample}")
        records.append({
            "sample": sample,
            "plate": plate,
            "cell_line_id": cell_line,
            "drug": drug,
            "drugname_drugconc": dose,
            "n_cells": int(state["n_cells"]),
            "n_cells_minimal": int(state["n_cells_minimal"]),
            "n_cells_obs_total": int(state["n_cells_obs_total"]),
            "mean_gene_count": state["gene_sum"] / state["gene_n"] if state["gene_n"] else source.get("mean_gene_count"),
            "mean_tscp_count": source.get("mean_tscp_count"),
            "mean_pcnt_mito": state["mito_sum"] / state["mito_n"] if state["mito_n"] else source.get("mean_pcnt_mito"),
        })
    return pd.DataFrame.from_records(records)


def stream_obs_metadata_sample_index(
    *,
    max_metadata_bytes: int,
    revision: str = DATASET_REVISION,
) -> pd.DataFrame:
    """Stream projected obs_metadata under a checked source-file byte ceiling."""
    try:
        from datasets import load_dataset
        from huggingface_hub import HfApi
    except ImportError as exc:
        raise RuntimeError("install the optional tahoe extra to scan obs_metadata") from exc
    if max_metadata_bytes <= 0:
        raise ValueError("max_metadata_bytes must be explicitly positive")
    info = HfApi().dataset_info(DATASET_ID, revision=revision, files_metadata=True)
    sibling = next(
        (item for item in (info.siblings or []) if item.rfilename == "metadata/obs_metadata.parquet"),
        None,
    )
    size = getattr(sibling, "size", None) if sibling is not None else None
    if size is None or size > max_metadata_bytes:
        raise ValueError(f"obs_metadata size {size} exceeds metadata byte budget {max_metadata_bytes}")
    sample_metadata = load_dataset(
        DATASET_ID, "sample_metadata", split="train", revision=revision,
    ).to_pandas()
    columns = [
        "sample", "plate", "drug", "drugname_drugconc", "cell_line",
        "gene_count", "pcnt_mito", "pass_filter",
    ]
    dataset = load_dataset(
        DATASET_ID, "obs_metadata", split="train", streaming=True, revision=revision,
    ).select_columns(columns)

    def projected_rows():
        for batch in dataset.iter(batch_size=50_000):
            for i in range(len(batch["sample"])):
                yield {column: batch[column][i] for column in columns}

    result = sample_index_from_metadata_rows(projected_rows(), sample_metadata)
    result.attrs["source_revision"] = getattr(info, "sha", revision)
    result.attrs["obs_metadata_source_bytes"] = int(size)
    result.attrs["metadata_scan_byte_budget"] = int(max_metadata_bytes)
    return result
""Pinned, metadata-only Hugging Face inspection helpers for Tahoe-100M."""
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
