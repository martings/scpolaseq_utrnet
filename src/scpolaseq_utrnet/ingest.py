"""PBMC ingestion with exact cell and PAS keys."""
from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution

from .config import Config
from .contracts import read_tsv, validate_cells, validate_counts


def _base(value: str) -> str:
    return value.split(".", 1)[0]


def _gene_id(config: Config, adata: ad.AnnData) -> tuple[str, str]:
    if config.gene.startswith("ENSG"):
        target = _base(config.gene)
        matches = adata.var[adata.var.gene_id_base.astype(str) == target]
    else:
        matches = adata.var[adata.var.gene_name.astype(str) == config.gene]
    ids = set(matches.gene_id_base.astype(str))
    if len(ids) != 1:
        raise ValueError(f"gene {config.gene!r} maps to {len(ids)} distinct gene IDs")
    return ids.pop(), str(matches.gene_name.iloc[0])


def load_cells(config: Config) -> tuple[pd.DataFrame, ad.AnnData]:
    annotations = read_tsv(config.upstream / "labels/cell_annotations.tsv")
    validate_cells(annotations)
    if set(annotations.sample_id) != {config.sample_id}:
        raise ValueError("sample ID mismatch")
    adata = ad.read_h5ad(config.curated_h5ad, backed="r")
    if not adata.obs_names.is_unique or "cell_id" not in adata.obs:
        raise ValueError("curated AnnData needs unique, explicit cell_id")
    obs = adata.obs.reset_index(drop=True)[["cell_id", "celltype_corrected", "qc_pass", "cluster_id"]].copy()
    if obs.cell_id.duplicated().any():
        raise ValueError("duplicate curated cell_id")
    cells = annotations.merge(obs, on="cell_id", how="left", validate="one_to_one", suffixes=("", "_curated"))
    cells["curated"] = cells.celltype_corrected.notna()
    if cells.loc[cells.curated, "qc_pass"].eq(False).any():
        raise ValueError("curated AnnData contains QC-fail cells")
    if not set(obs.cell_id).issubset(set(annotations.cell_id)):
        raise ValueError("curated cell absent from upstream annotations")
    return cells, adata


def load_pas_sites(config: Config, gene_id: str) -> pd.DataFrame:
    path = config.upstream / "pas_reference/pas_reference.tsv"
    parts = []
    for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=100_000):
        selected = chunk[chunk.gene_id.map(_base) == gene_id]
        if not selected.empty:
            parts.append(selected)
    if not parts:
        raise ValueError(f"no PAS reference sites for {gene_id}")
    sites = pd.concat(parts, ignore_index=True)
    if sites.site_id.duplicated().any() or sites.pas_reference_id.duplicated().any():
        raise ValueError("duplicate PAS reference identifier")
    sites["start"] = pd.to_numeric(sites.start, errors="raise").astype(int)
    sites["end"] = pd.to_numeric(sites.end, errors="raise").astype(int)
    if (sites.start != sites.end).any():
        raise ValueError("PAS reference is expected to contain 1-based cleavage points")
    if (~sites.strand.isin(["+", "-"])).any():
        raise ValueError("invalid PAS strand")
    return sites


def load_pas_counts(config: Config, sites: pd.DataFrame, cells: pd.DataFrame) -> pd.DataFrame:
    paths = sorted((config.upstream / "sierra_quant/cluster").glob("*.sierra_quant.tsv"))
    if not paths:
        raise FileNotFoundError("no cluster Sierra files")
    site_ids = set(sites.site_id)
    parts = []
    for path in paths:
        for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=200_000):
            selected = chunk[chunk.site_id.isin(site_ids)].copy()
            if not selected.empty:
                selected["source_file"] = str(path)
                parts.append(selected)
    counts = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=[
        "library_id", "group_level", "group_id", "site_id", "cell_barcode", "umi_count", "source_file"])
    validate_counts(counts, cells, sites)
    if not counts.empty and set(counts.group_level) != {"cluster"}:
        raise ValueError("mixed Sierra group levels in canonical observation view")
    key = cells[["sample_id", "library_id", "barcode_corrected", "cell_id", "curated", "celltype_corrected", "cluster_id"]]
    counts = counts.merge(key, left_on=["library_id", "cell_barcode"],
                          right_on=["library_id", "barcode_corrected"], how="left", validate="many_to_one")
    if counts.cell_id.isna().any() or (counts.group_id != counts.cluster_id).any():
        raise ValueError("Sierra cluster group disagrees with annotation cell cluster")
    counts["umi_count"] = pd.to_numeric(counts.umi_count).astype(int)
    return counts


def load_gene_expression(adata: ad.AnnData, gene_id: str) -> pd.DataFrame:
    var = adata.var
    index = np.flatnonzero(var.gene_id_base.astype(str).to_numpy() == gene_id)
    if len(index) != 1:
        raise ValueError("gene expression is absent or ambiguously mapped")
    if "counts" not in adata.layers:
        raise ValueError("curated AnnData lacks raw counts layer")
    column = adata.layers["counts"][:, int(index[0])]
    if hasattr(column, "toarray"):
        column = column.toarray()
    values = np.asarray(column).reshape(-1)
    return pd.DataFrame({"cell_id": adata.obs.cell_id.astype(str).to_numpy(),
                         "gene_id": gene_id, "value": values,
                         "layer": "counts", "normalization": "raw_UMI"})


def make_profiles(counts: pd.DataFrame, cells: pd.DataFrame, sites: pd.DataFrame,
                  gene_id: str) -> pd.DataFrame:
    curated = cells[cells.curated].copy()
    sizes = curated.groupby("celltype_corrected", observed=True).size().rename("n_cells")
    frame = counts[counts.curated == True].groupby(["celltype_corrected", "site_id"], as_index=False, observed=True).agg(
        umi_sum=("umi_count", "sum"), n_cells_with_umi=("cell_id", "nunique"))
    idx = pd.MultiIndex.from_product([sizes.index, sites.site_id], names=["celltype_corrected", "site_id"])
    frame = frame.set_index(["celltype_corrected", "site_id"]).reindex(idx, fill_value=0).reset_index()
    frame = frame.merge(sizes.reset_index(), on="celltype_corrected", validate="many_to_one")
    frame["total_gene_umi"] = frame.groupby("celltype_corrected", observed=True).umi_sum.transform("sum")
    frame["pas_proportion"] = np.where(frame.total_gene_umi > 0,
        frame.umi_sum / frame.total_gene_umi, np.nan)
    # Jeffreys Dirichlet marginal: UMI sampling interval, not biological replication.
    n_sites = len(sites)
    alpha = frame.umi_sum.astype(float) + 0.5
    beta = frame.total_gene_umi.astype(float) - frame.umi_sum.astype(float) + 0.5 * max(n_sites - 1, 1)
    lower = beta_distribution.ppf(0.025, alpha, beta)
    upper = beta_distribution.ppf(0.975, alpha, beta)
    frame["usage_ci_low"] = np.where(frame.total_gene_umi > 0,
                                      np.minimum(lower, frame.pas_proportion), np.nan)
    frame["usage_ci_high"] = np.where(frame.total_gene_umi > 0,
                                       np.maximum(upper, frame.pas_proportion), np.nan)
    frame["gene_id"] = gene_id
    frame["proportion_method"] = "pooled_curated_cell_UMI"
    frame["coverage_status"] = np.where(frame.total_gene_umi > 0, "observed", "no_gene_umi")
    return frame


def ingest_apafield_metrics(config: Config, gene_id: str) -> pd.DataFrame:
    """Ingest physical Apafield metrics (SI, magnetization, entropy) from downstream tables."""
    cols = [
        "celltype", "si", "magnetization_rank", "magnetization_genomic",
        "entropy", "entropy_normalized", "polarization_rank", "polarization_genomic"
    ]
    si_path = config.si_magnetization_table
    if not si_path.exists():
        return pd.DataFrame(columns=cols)

    base_id = _base(gene_id)
    df = pd.read_csv(si_path, sep="\t")
    filtered = df[df["gene_id_base"].astype(str).map(_base) == base_id].copy()
    if filtered.empty and "gene_name" in df.columns:
        filtered = df[df["gene_name"].astype(str).str.upper() == config.gene.upper()].copy()

    if filtered.empty:
        return pd.DataFrame(columns=cols)

    present = [c for c in cols if c in filtered.columns]
    return filtered[present].reset_index(drop=True)


def ingest_apa_priority(config: Config, gene_id: str) -> pd.DataFrame:
    """Ingest differential APA comparisons and priority rankings from downstream tables."""
    cols = [
        "comparison", "group_a", "group_b", "celltype", "delta_si",
        "delta_m_genomic_fixed", "adj_p_value", "is_robust",
        "priority_score", "endpoint_artifact_suspect"
    ]
    prio_path = config.apa_priority_table
    if not prio_path.exists():
        return pd.DataFrame(columns=cols)

    base_id = _base(gene_id)
    df = pd.read_csv(prio_path, sep="\t")
    filtered = df[df["gene_id_base"].astype(str).map(_base) == base_id].copy()
    if filtered.empty and "gene_name" in df.columns:
        filtered = df[df["gene_name"].astype(str).str.upper() == config.gene.upper()].copy()

    if filtered.empty:
        return pd.DataFrame(columns=cols)

    present = [c for c in cols if c in filtered.columns]
    return filtered[present].reset_index(drop=True)


def ingest_gene(config: Config) -> dict[str, pd.DataFrame]:
    cells, adata = load_cells(config)
    try:
        gene_id, gene_name = _gene_id(config, adata)
        sites = load_pas_sites(config, gene_id)
        counts = load_pas_counts(config, sites, cells)
        expression = load_gene_expression(adata, gene_id)
    finally:
        adata.file.close()
    profiles = make_profiles(counts, cells, sites, gene_id)
    sites["gene_name"] = gene_name

    apafield = ingest_apafield_metrics(config, gene_id)
    if not apafield.empty and "celltype" in apafield.columns:
        merge_cols = [c for c in apafield.columns if c != "celltype"]
        profiles = profiles.merge(
            apafield, left_on="celltype_corrected", right_on="celltype", how="left"
        )
        if "celltype" in profiles.columns:
            profiles = profiles.drop(columns=["celltype"])

    priority = ingest_apa_priority(config, gene_id)

    return {
        "cells": cells,
        "pas_sites": sites,
        "pas_counts": counts,
        "gene_expression": expression,
        "apa_profiles": profiles,
        "apafield_metrics": apafield,
        "apa_priority": priority,
    }


def write_ingestion(tables: dict[str, pd.DataFrame], directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for name, frame in tables.items():
        frame.to_csv(directory / f"{name}.csv", index=False)


def audit_gene_counts(config: Config, tables: dict[str, pd.DataFrame]) -> dict:
    """Compare the two Sierra grouping views at exact library/barcode/PAS keys."""
    target = set(tables["pas_sites"].site_id)
    parts = []
    files = sorted((config.upstream / "sierra_quant/cell_type").glob("*.sierra_quant.tsv"))
    if not files:
        return {"status": "missing_cell_type_grouping"}
    for path in files:
        for chunk in pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, chunksize=200_000):
            selected = chunk[chunk.site_id.isin(target)]
            if not selected.empty:
                parts.append(selected)
    other = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=tables["pas_counts"].columns)
    key = ["library_id", "cell_barcode", "site_id"]
    canonical = tables["pas_counts"].groupby(key).umi_count.sum().sort_index()
    other["umi_count"] = pd.to_numeric(other.umi_count, errors="raise").astype(int)
    grouped = other.groupby(key).umi_count.sum().sort_index()
    identical = canonical.equals(grouped)
    return {"status": "pass" if identical else "mismatch", "exact_key_equal": identical,
            "cluster_rows": len(tables["pas_counts"]), "cell_type_rows": len(other),
            "cluster_umi": int(tables["pas_counts"].umi_count.sum()), "cell_type_umi": int(other.umi_count.sum()),
            "distinct_observation_keys": len(canonical)}
