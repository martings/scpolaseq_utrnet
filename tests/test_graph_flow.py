"""Unit tests for discrete graph flows, node divergence, and trajectory contrasts."""
import pytest
import pandas as pd
import numpy as np

from scpolaseq_utrnet.flow import (
    compute_edge_flows,
    compute_node_divergence,
    compare_flow_divergence_trajectories,
)


@pytest.fixture
def synthetic_scores():
    # 2 target genes, 2 regulators, across 4 states (DF, DF_wk2, DHF, DHF_wk2) in Monocytes
    return pd.DataFrame([
        # State DF (day1)
        {"gene_name": "GENE1", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF", "interaction_score": 1.0},
        {"gene_name": "GENE2", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF", "interaction_score": 0.5},
        {"gene_name": "GENE1", "regulator": "TTP", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF", "interaction_score": 0.8},

        # State DF_wk2
        {"gene_name": "GENE1", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF_wk2", "interaction_score": 1.2},
        {"gene_name": "GENE2", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF_wk2", "interaction_score": 0.6},
        {"gene_name": "GENE1", "regulator": "TTP", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DF_wk2", "interaction_score": 0.7},

        # State DHF (day1)
        {"gene_name": "GENE1", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF", "interaction_score": 0.4},
        {"gene_name": "GENE2", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF", "interaction_score": 0.2},
        {"gene_name": "GENE1", "regulator": "TTP", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF", "interaction_score": 0.3},

        # State DHF_wk2
        {"gene_name": "GENE1", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF_wk2", "interaction_score": 1.5},
        {"gene_name": "GENE2", "regulator": "HuR", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF_wk2", "interaction_score": 0.8},
        {"gene_name": "GENE1", "regulator": "TTP", "regulator_class": "ARE-BP", "celltype_corrected": "Mono", "sample_id": "DHF_wk2", "interaction_score": 0.9},
    ])


def test_edge_flows_and_divergence(synthetic_scores):
    # DF trajectory: DF -> DF_wk2
    df_flows = compute_edge_flows(
        synthetic_scores, start_state="DF", end_state="DF_wk2", state_col="sample_id"
    )
    assert len(df_flows) == 3

    # HuR on GENE1 in DF: 1.2 - 1.0 = +0.2
    hur_g1 = df_flows[(df_flows["regulator"] == "HuR") & (df_flows["gene_name"] == "GENE1")].iloc[0]
    assert np.isclose(hur_g1["edge_flow"], 0.2)

    # Compute node divergence in DF
    div_df = compute_node_divergence(df_flows)
    hur_div = div_df[(div_df["node_id"] == "HuR") & (div_df["node_type"] == "regulator")].iloc[0]
    # Total HuR divergence = (1.2-1.0) + (0.6-0.5) = +0.3
    assert np.isclose(hur_div["field_divergence"], 0.3)
    assert hur_div["node_field_role"] == "sink_sponge_accumulation"


def test_trajectory_divergence_comparison(synthetic_scores):
    df_flows = compute_edge_flows(
        synthetic_scores, start_state="DF", end_state="DF_wk2", state_col="sample_id"
    )
    dhf_flows = compute_edge_flows(
        synthetic_scores, start_state="DHF", end_state="DHF_wk2", state_col="sample_id"
    )

    all_flows = pd.concat([df_flows, dhf_flows], ignore_index=True)
    divergence_df = compute_node_divergence(all_flows)

    contrast = compare_flow_divergence_trajectories(
        divergence_df, mild_trajectory=("DF", "DF_wk2"), severe_trajectory=("DHF", "DHF_wk2")
    )

    hur_contrast = contrast[(contrast["node_id"] == "HuR") & (contrast["node_type"] == "regulator")].iloc[0]
    # DF HuR div: 0.3
    # DHF HuR div: (1.5-0.4) + (0.8-0.2) = 1.1 + 0.6 = 1.7
    # Delta divergence: 1.7 - 0.3 = +1.4
    assert np.isclose(hur_contrast["div_mild"], 0.3)
    assert np.isclose(hur_contrast["div_severe"], 1.7)
    assert np.isclose(hur_contrast["delta_divergence"], 1.4)
    assert hur_contrast["divergence_phenotype"] == "severe_hyperaccumulation"
