"""Graph Field Theory: Discrete edge flows, node divergence, and regulatory flux modeling."""
from __future__ import annotations

import logging
from pathlib import Path
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_edge_flows(
    scores_df: pd.DataFrame,
    start_state: str,
    end_state: str,
    state_col: str = "sample_id",
    delta_t: float = 1.0,
) -> pd.DataFrame:
    """Compute discrete regulatory edge flows J_e = (W_e(end) - W_e(start)) / delta_t.

    Analogous to 1-form gauge field flows in lattice field theory.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Interaction scores table with columns:
        ['gene_name', 'regulator', 'regulator_class', 'celltype_corrected',
         state_col, 'interaction_score'].
    start_state : str
        Initial observation state (e.g., 'DF' or 'DHF').
    end_state : str
        Final observation state (e.g., 'DF_wk2' or 'DHF_wk2').
    state_col : str
        Column identifying the temporal state.
    delta_t : float
        Elapsed time interval.

    Returns
    -------
    pd.DataFrame
        Edge flow table:
        ['gene_name', 'regulator', 'regulator_class', 'celltype_corrected',
         'start_state', 'end_state', 'weight_start', 'weight_end',
         'edge_flow', 'flow_density']
    """
    if scores_df.empty:
        return pd.DataFrame(columns=[
            "gene_name", "regulator", "regulator_class", "celltype_corrected",
            "start_state", "end_state", "weight_start", "weight_end",
            "edge_flow", "flow_density"
        ])

    agg = scores_df.groupby(
        ["gene_name", "regulator", "regulator_class", "celltype_corrected", state_col],
        as_index=False
    )["interaction_score"].sum()

    sub_start = agg[agg[state_col] == start_state].set_index(
        ["gene_name", "regulator", "regulator_class", "celltype_corrected"]
    )
    sub_end = agg[agg[state_col] == end_state].set_index(
        ["gene_name", "regulator", "regulator_class", "celltype_corrected"]
    )

    all_keys = sub_start.index.union(sub_end.index)
    records = []

    for key in all_keys:
        gene_name, reg, reg_cls, ct = key
        w0 = float(sub_start.loc[key, "interaction_score"]) if key in sub_start.index else 0.0
        w1 = float(sub_end.loc[key, "interaction_score"]) if key in sub_end.index else 0.0
        j_e = (w1 - w0) / delta_t if delta_t > 0 else (w1 - w0)

        # Flow density relative to mean intensity
        mean_w = (w0 + w1) / 2.0
        f_dens = j_e / (mean_w + 1e-4)

        records.append({
            "gene_name": gene_name,
            "regulator": reg,
            "regulator_class": reg_cls,
            "celltype_corrected": ct,
            "start_state": start_state,
            "end_state": end_state,
            "weight_start": w0,
            "weight_end": w1,
            "edge_flow": j_e,
            "flow_density": f_dens,
        })

    return pd.DataFrame(records)


def compute_node_divergence(
    edge_flows: pd.DataFrame,
) -> pd.DataFrame:
    """Compute discrete field divergence (nabla . J)_v across regulator and gene nodes.

    For regulator r: (nabla . J)_r = sum_{g in N(r)} J_{(g, r)}
    Measures net accumulation (sink/sponge) or net depletion (source/release).

    Parameters
    ----------
    edge_flows : pd.DataFrame
        Output of `compute_edge_flows`.

    Returns
    -------
    pd.DataFrame
        Node divergence table:
        ['node_id', 'node_type', 'regulator_class', 'celltype_corrected',
         'start_state', 'end_state', 'field_divergence', 'abs_flux_sum',
         'node_field_role']
    """
    if edge_flows.empty:
        return pd.DataFrame(columns=[
            "node_id", "node_type", "regulator_class", "celltype_corrected",
            "start_state", "end_state", "field_divergence", "abs_flux_sum",
            "node_field_role"
        ])

    rows = []
    # 1. Regulator node divergence
    reg_grp = edge_flows.groupby(
        ["regulator", "regulator_class", "celltype_corrected", "start_state", "end_state"]
    )
    for (reg, reg_cls, ct, s0, s1), grp in reg_grp:
        div = float(grp["edge_flow"].sum())
        abs_sum = float(grp["edge_flow"].abs().sum())

        if div > 0.05:
            role = "sink_sponge_accumulation"
        elif div < -0.05:
            role = "source_regulator_liberation"
        else:
            role = "conservative_neutral_flux"

        rows.append({
            "node_id": reg,
            "node_type": "regulator",
            "regulator_class": reg_cls,
            "celltype_corrected": ct,
            "start_state": s0,
            "end_state": s1,
            "field_divergence": div,
            "abs_flux_sum": abs_sum,
            "node_field_role": role,
        })

    # 2. Gene node divergence
    gene_grp = edge_flows.groupby(
        ["gene_name", "celltype_corrected", "start_state", "end_state"]
    )
    for (gene, ct, s0, s1), grp in gene_grp:
        div = float(grp["edge_flow"].sum())
        abs_sum = float(grp["edge_flow"].abs().sum())

        if div > 0.05:
            role = "target_gain_of_regulation"
        elif div < -0.05:
            role = "target_shedding_of_regulation"
        else:
            role = "conservative_neutral_flux"

        rows.append({
            "node_id": gene,
            "node_type": "target_gene",
            "regulator_class": "target",
            "celltype_corrected": ct,
            "start_state": s0,
            "end_state": s1,
            "field_divergence": div,
            "abs_flux_sum": abs_sum,
            "node_field_role": role,
        })

    return pd.DataFrame(rows)


def compare_flow_divergence_trajectories(
    divergence_df: pd.DataFrame,
    mild_trajectory: tuple[str, str] = ("DF", "DF_wk2"),
    severe_trajectory: tuple[str, str] = ("DHF", "DHF_wk2"),
) -> pd.DataFrame:
    """Compare discrete divergence between mild and severe clinical trajectories.

    Delta (nabla . J) = (nabla . J)_{severe} - (nabla . J)_{mild}
    Isolates pathological source/sink anomalies in severe disease.

    Parameters
    ----------
    divergence_df : pd.DataFrame
        Output of `compute_node_divergence`.
    mild_trajectory : tuple[str, str]
        (start, end) for mild disease (e.g. ('DF', 'DF_wk2')).
    severe_trajectory : tuple[str, str]
        (start, end) for severe disease (e.g. ('DHF', 'DHF_wk2')).

    Returns
    -------
    pd.DataFrame
        Trajectory divergence contrast table:
        ['node_id', 'node_type', 'regulator_class', 'celltype_corrected',
         'div_mild', 'div_severe', 'delta_divergence', 'divergence_phenotype']
    """
    if divergence_df.empty:
        return pd.DataFrame(columns=[
            "node_id", "node_type", "regulator_class", "celltype_corrected",
            "div_mild", "div_severe", "delta_divergence", "divergence_phenotype"
        ])

    m0, m1 = mild_trajectory
    s0, s1 = severe_trajectory

    sub_mild = divergence_df[
        (divergence_df["start_state"] == m0) & (divergence_df["end_state"] == m1)
    ].set_index(["node_id", "node_type", "regulator_class", "celltype_corrected"])

    sub_sev = divergence_df[
        (divergence_df["start_state"] == s0) & (divergence_df["end_state"] == s1)
    ].set_index(["node_id", "node_type", "regulator_class", "celltype_corrected"])

    all_keys = sub_mild.index.union(sub_sev.index)
    records = []

    for key in all_keys:
        nid, ntype, reg_cls, ct = key
        d_mild = float(sub_mild.loc[key, "field_divergence"]) if key in sub_mild.index else 0.0
        d_sev = float(sub_sev.loc[key, "field_divergence"]) if key in sub_sev.index else 0.0
        delta_div = d_sev - d_mild

        if delta_div > 0.1:
            pheno = "severe_hyperaccumulation"
        elif delta_div < -0.1:
            pheno = "severe_depletion_collapse"
        else:
            pheno = "concordant_divergence"

        records.append({
            "node_id": nid,
            "node_type": ntype,
            "regulator_class": reg_cls,
            "celltype_corrected": ct,
            "div_mild": d_mild,
            "div_severe": d_sev,
            "delta_divergence": delta_div,
            "divergence_phenotype": pheno,
        })

    return pd.DataFrame(records)
