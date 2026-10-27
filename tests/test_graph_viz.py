"""Unit tests for cell-type regulatory network visualization and contrast graphs."""
import pytest
import pandas as pd
from pathlib import Path

from scpolaseq_utrnet.graph_viz import (
    filter_top_regulatory_edges,
    compute_contrast_edge_set,
    render_celltype_contrast_network,
    export_interactive_celltype_html,
)


@pytest.fixture
def mock_network_scores():
    return pd.DataFrame([
        {"gene_name": "ANXA1", "regulator": "HuR", "regulator_class": "ARE-BP", "population": "Mono", "interaction_score": 1.2},
        {"gene_name": "ANXA1", "regulator": "TTP", "regulator_class": "ARE-BP", "population": "Mono", "interaction_score": 0.8},
        {"gene_name": "PPIA", "regulator": "PTBP1", "regulator_class": "RRM", "population": "Mono", "interaction_score": 0.9},
        {"gene_name": "ANXA1", "regulator": "HuR", "regulator_class": "ARE-BP", "population": "Bcell", "interaction_score": 0.2},
        {"gene_name": "ANXA1", "regulator": "TTP", "regulator_class": "ARE-BP", "population": "Bcell", "interaction_score": 0.1},
        {"gene_name": "PPIA", "regulator": "PTBP1", "regulator_class": "RRM", "population": "Bcell", "interaction_score": 1.5},
    ])


def test_filter_top_regulatory_edges(mock_network_scores):
    edges = filter_top_regulatory_edges(mock_network_scores, population="Mono", top_n_edges=2)
    assert len(edges) == 2
    assert edges.iloc[0]["regulator"] == "HuR"
    assert edges.iloc[0]["interaction_score"] == 1.2


def test_compute_contrast_edge_set(mock_network_scores):
    edges_a, edges_b, edges_diff = compute_contrast_edge_set(
        mock_network_scores, pop_a="Mono", pop_b="Bcell", top_n_edges=5
    )
    assert not edges_a.empty
    assert not edges_b.empty
    assert not edges_diff.empty

    # HuR on ANXA1 in Bcell is 0.2, in Mono is 1.2 -> delta is -1.0
    hur_diff = edges_diff[(edges_diff["regulator"] == "HuR") & (edges_diff["gene_name"] == "ANXA1")].iloc[0]
    assert pytest.approx(hur_diff["delta_W"], 0.01) == -1.0


def test_render_contrast_network_file(mock_network_scores, tmp_path):
    out_png = tmp_path / "test_contrast.png"
    render_celltype_contrast_network(
        mock_network_scores, pop_a="Mono", pop_b="Bcell",
        output_path=out_png, label_a="Monocytes", label_b="B Cells"
    )
    assert out_png.exists()
    assert out_png.stat().st_size > 5000


def test_export_interactive_html(mock_network_scores, tmp_path):
    out_html = tmp_path / "test_interactive.html"
    export_interactive_celltype_html(
        mock_network_scores, pop_a="Mono", pop_b="Bcell",
        output_html=out_html, label_a="Monocytes", label_b="B Cells"
    )
    assert out_html.exists()
    assert "Plotly" in out_html.read_text()
