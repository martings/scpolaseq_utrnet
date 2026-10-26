"""Unit and integration tests for Graph Rewiring Velocity and dynamic trajectory modeling."""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.dynamics import (
    compute_graph_velocity,
    trajectory_apafield_velocity,
    export_trajectory_analysis,
)
from scpolaseq_utrnet.pipeline import build_gene_graphs
from scpolaseq_utrnet.visualization import render_gene


def test_compute_graph_velocity():
    # 3 states in trajectory: State1 -> State2 -> State3
    # RegA: weight 0.1 -> 0.4 -> 0.8 (monotone_gain)
    # RegB: weight 0.8 -> 0.4 -> 0.1 (monotone_loss)
    # RegC: weight 0.1 -> 0.9 -> 0.2 (transient_activation)
    # RegD: weight 0.5 -> 0.5 -> 0.5 (stable)
    scores = pd.DataFrame([
        # State1
        {"population": "State1", "regulator": "RegA", "regulator_class": "RBP", "genomic_region": "reg1", "interaction_score": 0.1},
        {"population": "State1", "regulator": "RegB", "regulator_class": "RBP", "genomic_region": "reg2", "interaction_score": 0.8},
        {"population": "State1", "regulator": "RegC", "regulator_class": "RBP", "genomic_region": "reg3", "interaction_score": 0.1},
        {"population": "State1", "regulator": "RegD", "regulator_class": "RBP", "genomic_region": "reg4", "interaction_score": 0.5},
        # State2
        {"population": "State2", "regulator": "RegA", "regulator_class": "RBP", "genomic_region": "reg1", "interaction_score": 0.4},
        {"population": "State2", "regulator": "RegB", "regulator_class": "RBP", "genomic_region": "reg2", "interaction_score": 0.4},
        {"population": "State2", "regulator": "RegC", "regulator_class": "RBP", "genomic_region": "reg3", "interaction_score": 0.9},
        {"population": "State2", "regulator": "RegD", "regulator_class": "RBP", "genomic_region": "reg4", "interaction_score": 0.5},
        # State3
        {"population": "State3", "regulator": "RegA", "regulator_class": "RBP", "genomic_region": "reg1", "interaction_score": 0.8},
        {"population": "State3", "regulator": "RegB", "regulator_class": "RBP", "genomic_region": "reg2", "interaction_score": 0.1},
        {"population": "State3", "regulator": "RegC", "regulator_class": "RBP", "genomic_region": "reg3", "interaction_score": 0.2},
        {"population": "State3", "regulator": "RegD", "regulator_class": "RBP", "genomic_region": "reg4", "interaction_score": 0.5},
    ])

    order = ["State1", "State2", "State3"]
    vel = compute_graph_velocity(scores, order)
    assert not vel.empty
    assert len(vel) == 8  # 4 regulators * 2 steps

    # Check RegA (monotone gain)
    rega_step1 = vel[(vel.regulator == "RegA") & (vel.step == 1)].iloc[0]
    assert rega_step1.velocity == pytest.approx(0.3)
    assert rega_step1.trajectory_pattern == "monotone_gain"

    # Check RegB (monotone loss)
    regb_step1 = vel[(vel.regulator == "RegB") & (vel.step == 1)].iloc[0]
    assert regb_step1.velocity == pytest.approx(-0.4)
    assert regb_step1.trajectory_pattern == "monotone_loss"

    # Check RegC (transient activation)
    regc = vel[vel.regulator == "RegC"].iloc[0]
    assert regc.trajectory_pattern == "transient_activation"

    # Check RegD (stable)
    regd = vel[vel.regulator == "RegD"].iloc[0]
    assert regd.trajectory_pattern == "stable"


def test_trajectory_apafield_velocity():
    apa_df = pd.DataFrame([
        {"celltype": "S1", "si": 0.2, "magnetization_genomic": -0.6},
        {"celltype": "S2", "si": 0.5, "magnetization_genomic": 0.0},
        {"celltype": "S3", "si": 0.9, "magnetization_genomic": 0.7},
    ])
    vel = trajectory_apafield_velocity(apa_df, ["S1", "S2", "S3"])
    assert len(vel) == 2
    assert vel.iloc[0].delta_si == pytest.approx(0.3)
    assert vel.iloc[1].delta_si == pytest.approx(0.4)
    assert vel.iloc[0].delta_mag == pytest.approx(0.6)


def test_export_and_render_end_to_end(tmp_path: Path):
    config_path = make_fixture(tmp_path)
    config = load_config(config_path)

    result = build_gene_graphs(config)
    directory = result["directory"]

    # Test trajectory export
    traj_res = export_trajectory_analysis(directory, ["A", "B"])
    assert traj_res["velocity"].exists()

    # Test render_gene
    html, static = render_gene(directory)
    assert html.exists() and html.stat().st_size > 1000
    assert static.exists() and static.stat().st_size > 1000
