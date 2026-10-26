"""Measured mRNA regulator abundance; miRNA is prior-only without small-RNA data."""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from .config import Config


def regulator_abundance(config: Config, regulatory_sites: pd.DataFrame,
                        cells: pd.DataFrame) -> pd.DataFrame:
    populations = cells.loc[cells.curated, ["cell_id", "celltype_corrected"]]
    rows = []
    if regulatory_sites.empty:
        return pd.DataFrame(columns=["regulator", "regulator_class", "population", "mean_umi_per_cell", "abundance_source", "uncertainty"])
    adata = ad.read_h5ad(config.curated_h5ad, backed="r")
    try:
        var = adata.var
        for regulator, kind in regulatory_sites[["regulator", "regulator_class"]].drop_duplicates().itertuples(index=False, name=None):
            index = np.flatnonzero(var.gene_name.astype(str).to_numpy() == regulator)
            measured = kind != "miRNA" and len(index) == 1 and "counts" in adata.layers
            if measured:
                column = adata.layers["counts"][:, int(index[0])]
                if hasattr(column, "toarray"):
                    column = column.toarray()
                values = np.asarray(column).reshape(-1)
                per_cell = pd.DataFrame({"cell_id": adata.obs.cell_id.astype(str).to_numpy(), "umi": values})
                means = populations.merge(per_cell, on="cell_id", validate="one_to_one").groupby("celltype_corrected", observed=True).umi.mean()
            else:
                means = {}
            for population in sorted(populations.celltype_corrected.unique()):
                rows.append({"regulator": regulator, "regulator_class": kind, "population": population,
                             "mean_umi_per_cell": float(means[population]) if measured else np.nan,
                             "abundance_source": "curated_10x_mRNA_counts" if measured else "missing_prior_only",
                             "uncertainty": "single_sample_cell_mean" if measured else "unmeasured"})
    finally:
        adata.file.close()
    return pd.DataFrame(rows)
