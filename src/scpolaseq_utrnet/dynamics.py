"""Graph Rewiring Velocity and dynamic trajectory modeling."""
from __future__ import annotations

import logging
from pathlib import Path
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def compute_graph_velocity(
    scores: pd.DataFrame,
    trajectory_order: list[str],
) -> pd.DataFrame:
    """Compute step-wise Graph Rewiring Velocity along an ordered trajectory of cell states.

    Velocity V_e(t) = (W_e(t+1) - W_e(t)) / delta_t
    Quantifies the directional rate of change in candidate regulatory interactions.
    """
    if scores.empty or len(trajectory_order) < 2:
        return pd.DataFrame(columns=[
            "regulator", "regulator_class", "genomic_region", "step",
            "source_state", "target_state", "weight_start", "weight_end",
            "velocity", "acceleration", "trajectory_pattern", "total_trajectory_delta"
        ])

    # Filter to states present in trajectory and scores
    available_states = [s for s in trajectory_order if s in scores.population.values]
    if len(available_states) < 2:
        logger.warning("Fewer than 2 trajectory states available in scores: %s", available_states)
        return pd.DataFrame(columns=[
            "regulator", "regulator_class", "genomic_region", "step",
            "source_state", "target_state", "weight_start", "weight_end",
            "velocity", "acceleration", "trajectory_pattern", "total_trajectory_delta"
        ])

    # Aggregate interaction score per state, regulator, region
    sub = scores[scores.population.isin(available_states)]
    agg = (
        sub.groupby(["population", "regulator", "regulator_class", "genomic_region"], as_index=False)
        ["interaction_score"]
        .sum()
    )

    # Pivot into state matrix: index = (regulator, regulator_class, genomic_region), columns = available_states
    matrix = agg.pivot(
        index=["regulator", "regulator_class", "genomic_region"],
        columns="population",
        values="interaction_score",
    ).reindex(columns=available_states).fillna(0.0)

    records = []
    n_steps = len(available_states) - 1

    for (reg, reg_cls, region), row in matrix.iterrows():
        weights = row.values.astype(float)
        velocities = np.diff(weights)
        total_delta = float(weights[-1] - weights[0])

        # Pattern classification
        if np.all(np.abs(velocities) < 1e-5):
            pattern = "stable"
        elif np.all(velocities >= -1e-5):
            pattern = "monotone_gain"
        elif np.all(velocities <= 1e-5):
            pattern = "monotone_loss"
        elif weights[1:-1].size > 0 and np.max(weights[1:-1]) > max(weights[0], weights[-1]) + 1e-4:
            pattern = "transient_activation"
        elif weights[1:-1].size > 0 and np.min(weights[1:-1]) < min(weights[0], weights[-1]) - 1e-4:
            pattern = "transient_repression"
        else:
            pattern = "dynamic_complex"

        accelerations = np.diff(velocities) if len(velocities) > 1 else np.array([0.0])

        for step_idx in range(n_steps):
            s_from = available_states[step_idx]
            s_to = available_states[step_idx + 1]
            v = float(velocities[step_idx])
            acc = float(accelerations[step_idx]) if step_idx < len(accelerations) else float("nan")

            records.append({
                "regulator": reg,
                "regulator_class": reg_cls,
                "genomic_region": region,
                "step": step_idx + 1,
                "source_state": s_from,
                "target_state": s_to,
                "weight_start": float(weights[step_idx]),
                "weight_end": float(weights[step_idx + 1]),
                "velocity": v,
                "acceleration": acc,
                "trajectory_pattern": pattern,
                "total_trajectory_delta": total_delta,
            })

    return pd.DataFrame(records)


def trajectory_apafield_velocity(
    apafield_df: pd.DataFrame,
    trajectory_order: list[str],
) -> pd.DataFrame:
    """Compute rate of change of physical APA metrics (SI, Magnetization) along a trajectory."""
    if apafield_df.empty or len(trajectory_order) < 2:
        return pd.DataFrame(columns=[
            "step", "source_state", "target_state", "si_start", "si_end",
            "delta_si", "velocity_si", "mag_start", "mag_end", "delta_mag", "velocity_mag"
        ])

    celltype_col = "celltype" if "celltype" in apafield_df.columns else "celltype_corrected"
    sub = apafield_df.set_index(celltype_col)

    available = [s for s in trajectory_order if s in sub.index]
    if len(available) < 2:
        return pd.DataFrame(columns=[
            "step", "source_state", "target_state", "si_start", "si_end",
            "delta_si", "velocity_si", "mag_start", "mag_end", "delta_mag", "velocity_mag"
        ])

    records = []
    for step_idx in range(len(available) - 1):
        s_from = available[step_idx]
        s_to = available[step_idx + 1]

        si_0 = float(sub.loc[s_from, "si"]) if "si" in sub.columns and pd.notna(sub.loc[s_from, "si"]) else np.nan
        si_1 = float(sub.loc[s_to, "si"]) if "si" in sub.columns and pd.notna(sub.loc[s_to, "si"]) else np.nan
        delta_si = si_1 - si_0 if pd.notna(si_0) and pd.notna(si_1) else np.nan

        mag_col = "magnetization_genomic" if "magnetization_genomic" in sub.columns else None
        mag_0 = float(sub.loc[s_from, mag_col]) if mag_col and pd.notna(sub.loc[s_from, mag_col]) else np.nan
        mag_1 = float(sub.loc[s_to, mag_col]) if mag_col and pd.notna(sub.loc[s_to, mag_col]) else np.nan
        delta_mag = mag_1 - mag_0 if pd.notna(mag_0) and pd.notna(mag_1) else np.nan

        records.append({
            "step": step_idx + 1,
            "source_state": s_from,
            "target_state": s_to,
            "si_start": si_0,
            "si_end": si_1,
            "delta_si": delta_si,
            "velocity_si": delta_si,  # delta_t = 1
            "mag_start": mag_0,
            "mag_end": mag_1,
            "delta_mag": delta_mag,
            "velocity_mag": delta_mag,
        })

    return pd.DataFrame(records)


def export_trajectory_analysis(
    directory: Path,
    trajectory_order: list[str],
) -> dict[str, Path]:
    """Execute trajectory velocity analysis and export tabular summaries."""
    scores_path = directory / "interaction_scores.csv"
    if not scores_path.exists():
        raise FileNotFoundError(f"interaction_scores.csv not found in {directory}")

    scores = pd.read_csv(scores_path)
    velocity_df = compute_graph_velocity(scores, trajectory_order)
    vel_path = directory / "trajectory_velocity.csv"
    velocity_df.to_csv(vel_path, index=False)

    apa_path = directory / "apafield_metrics.csv"
    apa_vel_path = directory / "trajectory_apafield.csv"
    if apa_path.exists():
        apa_df = pd.read_csv(apa_path)
        apa_vel = trajectory_apafield_velocity(apa_df, trajectory_order)
        apa_vel.to_csv(apa_vel_path, index=False)

    return {
        "velocity": vel_path,
        "apafield_velocity": apa_vel_path,
    }
