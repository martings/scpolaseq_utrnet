"""Unit and integration tests for ceRNA cross-talk and RBP-sponge network modeling."""
from __future__ import annotations

from pathlib import Path
import networkx as nx
import pandas as pd
import pytest

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.cerna import (
    compute_cerna_coupling,
    build_cerna_network,
    export_cerna_analysis,
)


def test_cerna_coupling_and_network_math():
    # 2 populations (PopA, PopB)
    # 3 genes: Gene1, Gene2, Gene3
    # 2 regulators: RegX, RegY
    # Gene1 and Gene2 both bind RegX (strong ceRNA pair)
    # Gene3 binds RegY only (uncoupled to Gene1/Gene2)
    effective_df = pd.DataFrame([
        # PopA
        {"population": "PopA", "gene": "Gene1", "regulator": "RegX", "regulator_class": "ARE-BP", "effective_sites": 2.0, "gene_expression": 10.0, "gene_load": 20.0},
        {"population": "PopA", "gene": "Gene2", "regulator": "RegX", "regulator_class": "ARE-BP", "effective_sites": 1.0, "gene_expression": 10.0, "gene_load": 10.0},
        {"population": "PopA", "gene": "Gene3", "regulator": "RegY", "regulator_class": "PUF", "effective_sites": 3.0, "gene_expression": 5.0, "gene_load": 15.0},
        # PopB
        {"population": "PopB", "gene": "Gene1", "regulator": "RegX", "regulator_class": "ARE-BP", "effective_sites": 1.0, "gene_expression": 5.0, "gene_load": 5.0},
    ])

    load_df = pd.DataFrame([
        {"population": "PopA", "regulator": "RegX", "regulator_class": "ARE-BP", "total_load": 30.0},
        {"population": "PopA", "regulator": "RegY", "regulator_class": "PUF", "total_load": 15.0},
        {"population": "PopB", "regulator": "RegX", "regulator_class": "ARE-BP", "total_load": 5.0},
    ])

    coupling = compute_cerna_coupling(effective_df, load_df, "PopA", min_coupling=0.01)
    assert not coupling.empty
    # Gene1 and Gene2 share RegX, so they should be coupled with score 1.0 (cosine similarity of parallel vectors)
    pair = coupling[(coupling.gene_a == "Gene1") & (coupling.gene_b == "Gene2")].iloc[0]
    assert pair.coupling_score == pytest.approx(1.0, rel=1e-3)
    assert pair.shared_regulators_count == 1
    assert pair.shared_regulators_list == "RegX"

    # Gene3 does not share regulators with Gene1 or Gene2, so no edges involving Gene3
    assert "Gene3" not in coupling.gene_a.values and "Gene3" not in coupling.gene_b.values

    # Build NetworkX graph
    graph = build_cerna_network(coupling, effective_df, load_df, "PopA")
    assert f"gene:Gene1" in graph
    assert f"gene:Gene2" in graph
    assert f"gene:Gene3" in graph
    assert f"regulator:RegX" in graph
    assert f"regulator:RegY" in graph

    # Titration edge between Gene1 and RegX
    assert graph.has_edge("gene:Gene1", "regulator:RegX")
    assert graph["gene:Gene1"]["regulator:RegX"]["kind"] == "titrates_regulator"

    # Crosstalk edge between Gene1 and Gene2
    assert graph.has_edge("gene:Gene1", "gene:Gene2")
    assert graph["gene:Gene1"]["gene:Gene2"]["kind"] == "cerna_crosstalk"


def test_cerna_export_end_to_end(tmp_path: Path):
    config_path = make_fixture(tmp_path)
    config = load_config(config_path)

    # Run ceRNA analysis on fixture gene PLUS
    res = export_cerna_analysis(config, genes=["PLUS"], output_dir=tmp_path / "cerna_out")
    assert res["effective_sites"].exists()
    assert res["cellular_load"].exists()
    assert res["coupling_matrix"].exists()
    assert res["graphs"].exists()
