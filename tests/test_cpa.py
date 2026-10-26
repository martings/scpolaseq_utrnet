"""Unit and integration tests for Cleavage & Polyadenylation (CPA) machinery modeling."""
from __future__ import annotations

from pathlib import Path
import anndata as ad
import networkx as nx
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.cpa import (
    CPA_MACHINERY_CATALOG,
    extract_cpa_expression,
    correlate_cpa_with_apa,
    build_cpa_regulatory_edges,
)
from scpolaseq_utrnet.graphs import build_graph
from scpolaseq_utrnet.pipeline import build_gene_graphs


def test_cpa_catalog_completeness():
    assert len(CPA_MACHINERY_CATALOG) >= 38
    symbols = {entry["symbol"] for entry in CPA_MACHINERY_CATALOG}
    required = {
        "NUDT21", "CPSF6", "CPSF7", "CPSF1", "CPSF2", "CPSF3", "CPSF4", "WDR33", "FIP1L1",
        "CSTF1", "CSTF2", "CSTF2T", "CSTF3", "PABPN1", "PCF11", "SRSF3", "PTBP1", "MATR3"
    }
    assert required.issubset(symbols)
    for entry in CPA_MACHINERY_CATALOG:
        assert entry["complex"] in {
            "CPSF", "CSTF", "CFIm", "CFIIm", "Polymerase", "PAB",
            "Cleavage_Scaffold", "Cleavage_Regulator", "Phosphatase", "Splicing_APA"
        }
        assert entry["pmid"] is not None


def test_extract_cpa_expression_mock():
    # Empty case
    empty_res = extract_cpa_expression(None)
    assert empty_res.empty

    # Create mock AnnData with 3 populations and 2 CPA factors
    obs = pd.DataFrame({
        "celltype_corrected": ["PopA", "PopA", "PopB", "PopB", "PopC", "PopC"],
        "qc_pass": [True, True, True, True, True, True],
    }, index=[f"cell_{i}" for i in range(6)])

    var = pd.DataFrame({
        "gene_id": ["ENSG_NUDT21.1", "ENSG_CSTF2.1"],
        "gene_name": ["NUDT21", "CSTF2"],
    }, index=["NUDT21", "CSTF2"])

    # PopA: NUDT21=10, CSTF2=0
    # PopB: NUDT21=5,  CSTF2=5
    # PopC: NUDT21=0,  CSTF2=10
    counts = np.array([
        [10, 0],
        [10, 0],
        [5, 5],
        [5, 5],
        [0, 10],
        [0, 10],
    ])
    adata = ad.AnnData(X=sp.csr_matrix(counts), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(counts)

    cpa_expr = extract_cpa_expression(adata, population_col="celltype_corrected")
    assert not cpa_expr.empty
    assert set(cpa_expr.population) == {"PopA", "PopB", "PopC"}
    assert set(cpa_expr.factor) == {"NUDT21", "CSTF2"}

    nudt21_a = cpa_expr[(cpa_expr.factor == "NUDT21") & (cpa_expr.population == "PopA")].iloc[0]
    assert nudt21_a.mean_expr == pytest.approx(10.0)
    assert nudt21_a.pct_expressing == pytest.approx(1.0)

    nudt21_c = cpa_expr[(cpa_expr.factor == "NUDT21") & (cpa_expr.population == "PopC")].iloc[0]
    assert nudt21_c.mean_expr == pytest.approx(0.0)
    assert nudt21_c.pct_expressing == pytest.approx(0.0)


def test_correlate_cpa_with_apa():
    # Build mock cpa_expr
    cpa_expr = pd.DataFrame([
        {"population": "PopA", "factor": "NUDT21", "mean_expr": 10.0, "complex": "CFIm", "prior_effect": "represses_proximal_pas", "pmid": "29107936"},
        {"population": "PopB", "factor": "NUDT21", "mean_expr": 5.0,  "complex": "CFIm", "prior_effect": "represses_proximal_pas", "pmid": "29107936"},
        {"population": "PopC", "factor": "NUDT21", "mean_expr": 0.0,  "complex": "CFIm", "prior_effect": "represses_proximal_pas", "pmid": "29107936"},
    ])

    # PopA is distal (high SI), PopC is proximal (low SI)
    apa_profiles = pd.DataFrame([
        {"celltype_corrected": "PopA", "site_id": "site_distal", "pas_proportion": 0.9, "si": 0.9, "magnetization_genomic": 0.8},
        {"celltype_corrected": "PopA", "site_id": "site_proximal", "pas_proportion": 0.1, "si": 0.9, "magnetization_genomic": 0.8},
        {"celltype_corrected": "PopB", "site_id": "site_distal", "pas_proportion": 0.5, "si": 0.5, "magnetization_genomic": 0.0},
        {"celltype_corrected": "PopB", "site_id": "site_proximal", "pas_proportion": 0.5, "si": 0.5, "magnetization_genomic": 0.0},
        {"celltype_corrected": "PopC", "site_id": "site_distal", "pas_proportion": 0.1, "si": 0.1, "magnetization_genomic": -0.8},
        {"celltype_corrected": "PopC", "site_id": "site_proximal", "pas_proportion": 0.9, "si": 0.1, "magnetization_genomic": -0.8},
    ])

    correlations = correlate_cpa_with_apa(cpa_expr, apa_profiles)
    assert not correlations.empty

    # NUDT21 vs SI should be perfectly positively correlated (r = 1.0)
    nudt_si = correlations[(correlations.factor == "NUDT21") & (correlations.target_type == "metric_si")].iloc[0]
    assert nudt_si.spearman_r == pytest.approx(1.0)
    assert nudt_si.concordant_with_prior is True

    # NUDT21 vs proximal site should be perfectly negatively correlated (r = -1.0)
    nudt_prox = correlations[(correlations.factor == "NUDT21") & (correlations.target_id == "site_proximal")].iloc[0]
    assert nudt_prox.spearman_r == pytest.approx(-1.0)

    # Build regulatory edges
    edges = build_cpa_regulatory_edges(correlations, r_thresh=0.5, p_thresh=0.1)
    assert not edges.empty
    prox_edge = edges[edges.target == "pas:site_proximal"].iloc[0]
    assert prox_edge.source == "cpa:NUDT21"
    assert prox_edge.direction == "represses_usage"
    assert prox_edge.correlation_r == pytest.approx(-1.0)


def test_cpa_graph_integration(tmp_path: Path):
    config_path = make_fixture(tmp_path)
    config = load_config(config_path)

    # Add NUDT21 to the fixture's adata_celltyped.h5ad
    adata = ad.read_h5ad(config.curated_h5ad)
    new_var = adata.var.copy()
    new_var.loc["NUDT21"] = ["ENSGNUDT21.1", "ENSGNUDT21", "NUDT21"]
    
    # 2 cells: cell1 has 5 counts, cell2 has 0 counts
    old_counts = adata.layers["counts"].toarray()
    extra_col = np.array([[5], [0]])
    new_counts = np.hstack([old_counts, extra_col])
    
    new_adata = ad.AnnData(X=sp.csr_matrix(new_counts), obs=adata.obs.copy(), var=new_var)
    new_adata.layers["counts"] = sp.csr_matrix(new_counts)
    new_adata.write_h5ad(config.curated_h5ad)

    # Run build_gene_graphs
    result = build_gene_graphs(config)
    directory = result["directory"]

    assert (directory / "cpa_expression.csv").exists()
    assert (directory / "cpa_correlations.csv").exists()
    assert (directory / "cpa_regulatory_edges.csv").exists()

    cpa_expr = pd.read_csv(directory / "cpa_expression.csv")
    assert "NUDT21" in cpa_expr.factor.values

    # Check graphml contains cpa node and edge if present
    graph_files = list((directory / "graphs").glob("*.graphml"))
    assert len(graph_files) > 0
    g = nx.read_graphml(graph_files[0])
    assert any("population" in str(n) for n in g.nodes)
