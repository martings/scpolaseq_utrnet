"""Explanatory, publication-grade regulatory network visualization by cell type and comparative contrast."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path as MplPath
import plotly.graph_objects as go

logger = logging.getLogger(__name__)

# Standardized color scheme for regulatory classes
CLASS_PALETTE = {
    "ARE-BP": "#7c3aed",        # Purple
    "KH_domain": "#0284c7",     # Cerulean
    "RRM": "#2563eb",           # Royal Blue
    "Zinc_knuckle": "#ea580c",   # Orange
    "Zinc_finger": "#ec4899",   # Pink
    "hnRNP": "#4338ca",         # Indigo
    "CSD": "#0d9488",           # Teal
    "Endonuclease": "#dc2626",  # Red
    "Helicase": "#d97706",      # Amber
    "Deadenylase": "#b45309",   # Brown
    "m6A_adaptor": "#9333ea",   # Violet
    "CPA_machinery": "#059669", # Emerald
    "miRNA": "#f43f5e",         # Rose
    "other": "#64748b",         # Slate
}


def filter_top_regulatory_edges(
    scores_df: pd.DataFrame,
    population: str,
    top_n_edges: int = 35,
    min_score: float = 0.05,
    pop_col: str = "population",
) -> pd.DataFrame:
    """Filter to the most informative, legible subset of regulatory interactions for a population."""
    sub = scores_df[scores_df[pop_col] == population].copy()
    if sub.empty:
        return pd.DataFrame()

    agg = sub.groupby(
        ["gene_name", "regulator", "regulator_class"], as_index=False
    )["interaction_score"].sum()

    agg = agg[agg["interaction_score"] >= min_score]
    agg = agg.sort_values("interaction_score", ascending=False).head(top_n_edges)
    return agg.reset_index(drop=True)


def compute_contrast_edge_set(
    scores_df: pd.DataFrame,
    pop_a: str,
    pop_b: str,
    top_n_edges: int = 30,
    pop_col: str = "population",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Extract top interactions in pop_a, pop_b, and their differential rewiring delta."""
    sub = scores_df[scores_df[pop_col].isin([pop_a, pop_b])].copy()
    if sub.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    agg = sub.groupby(
        ["gene_name", "regulator", "regulator_class", pop_col], as_index=False
    )["interaction_score"].sum()

    piv = agg.pivot(
        index=["gene_name", "regulator", "regulator_class"],
        columns=pop_col,
        values="interaction_score"
    ).fillna(0.0).reset_index()

    if pop_a not in piv.columns:
        piv[pop_a] = 0.0
    if pop_b not in piv.columns:
        piv[pop_b] = 0.0

    piv["delta_W"] = piv[pop_b] - piv[pop_a]
    piv["abs_delta"] = piv["delta_W"].abs()

    # Edges for A, B, and Differential
    edges_a = piv[piv[pop_a] > 0.05].sort_values(pop_a, ascending=False).head(top_n_edges)
    edges_b = piv[piv[pop_b] > 0.05].sort_values(pop_b, ascending=False).head(top_n_edges)
    edges_diff = piv.sort_values("abs_delta", ascending=False).head(top_n_edges)

    return edges_a, edges_b, edges_diff


def _draw_bipartite_network_panel(
    ax: plt.Axes,
    edges: pd.DataFrame,
    weight_col: str,
    title: str,
    is_contrast: bool = False,
    si_mapping: dict[str, float] | None = None,
):
    """Render a clean, explanatory bipartite graph panel without edge congestion."""
    ax.set_facecolor("#ffffff")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    if edges.empty:
        ax.text(0.5, 0.5, "No regulatory interactions detected\nat current threshold",
                ha="center", va="center", fontsize=12, color="#94a3b8")
        ax.set_title(title, fontsize=12, fontweight="bold", pad=12)
        return

    genes = sorted(edges["gene_name"].unique())
    regs = sorted(edges["regulator"].unique())

    n_g = len(genes)
    n_r = len(regs)

    y_genes = np.linspace(0.9, 0.1, n_g) if n_g > 1 else np.array([0.5])
    y_regs = np.linspace(0.9, 0.1, n_r) if n_r > 1 else np.array([0.5])

    pos_genes = {g: (0.15, y) for g, y in zip(genes, y_genes)}
    pos_regs = {r: (0.85, y) for r, y in zip(regs, y_regs)}

    max_w = edges[weight_col].abs().max() if edges[weight_col].abs().max() > 0 else 1.0

    # Draw smooth bezier curve edges
    for _, row in edges.iterrows():
        g = row["gene_name"]
        r = row["regulator"]
        w = row[weight_col]

        x0, y0 = pos_genes[g]
        x1, y1 = pos_regs[r]

        # Spline control points
        ctrl1 = (x0 + 0.35, y0)
        ctrl2 = (x1 - 0.35, y1)

        path_data = [
            (MplPath.MOVETO, (x0, y0)),
            (MplPath.CURVE4, ctrl1),
            (MplPath.CURVE4, ctrl2),
            (MplPath.CURVE4, (x1, y1)),
        ]
        path = MplPath([p[1] for p in path_data], [p[0] for p in path_data])

        if is_contrast:
            # Contrast encoding: gained = orange/red, lost = blue/cyan
            if w > 0.1:
                color = "#d95f02"  # Gained in Pop B
                alpha = min(0.85, 0.35 + 0.5 * (w / max_w))
                lw = 1.0 + 3.5 * (w / max_w)
            elif w < -0.1:
                color = "#2b83ba"  # Lost in Pop B (specific to A)
                alpha = min(0.85, 0.35 + 0.5 * (abs(w) / max_w))
                lw = 1.0 + 3.5 * (abs(w) / max_w)
            else:
                color = "#cbd5e1"
                alpha = 0.25
                lw = 0.8
        else:
            # Single population encoding: color by regulator class
            reg_cls = row["regulator_class"]
            color = CLASS_PALETTE.get(reg_cls, "#64748b")
            norm_w = min(1.0, max(0.0, w / max_w))
            alpha = min(0.9, 0.25 + 0.65 * norm_w)
            lw = 0.8 + 3.2 * norm_w

        patch = patches.PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=alpha, zorder=2)
        ax.add_patch(patch)

    # Draw Gene nodes (Rectangles on left)
    for g, (x, y) in pos_genes.items():
        # Node color based on SI or shortening if available
        if si_mapping and g in si_mapping:
            si_val = si_mapping[g]
            gene_color = "#fef08a" if si_val > 0.6 else "#fed7aa" if si_val > 0.4 else "#bae6fd"
        else:
            gene_color = "#e2e8f0"

        rect = patches.FancyBboxPatch(
            (x - 0.10, y - 0.025), 0.10, 0.05,
            boxstyle="round,pad=0.015",
            facecolor=gene_color, edgecolor="#334155", lw=1.2, zorder=4
        )
        ax.add_patch(rect)
        ax.text(x - 0.05, y, g, ha="center", va="center", fontsize=9, fontweight="bold", color="#0f172a", zorder=5)

    # Draw Regulator nodes (Circles/Octagons on right)
    reg_class_map = edges.groupby("regulator")["regulator_class"].first().to_dict()
    for r, (x, y) in pos_regs.items():
        r_cls = reg_class_map.get(r, "other")
        r_color = CLASS_PALETTE.get(r_cls, "#64748b")

        circle = patches.Circle((x, y), 0.022, facecolor=r_color, edgecolor="#0f172a", lw=1.0, zorder=4)
        ax.add_patch(circle)
        ax.text(x + 0.035, y, r, ha="left", va="center", fontsize=8.5, fontweight="semibold", color="#1e293b", zorder=5)

    # Column Headers
    ax.text(0.10, 0.97, "Target Genes", ha="center", va="center", fontsize=10, fontweight="bold", color="#475569")
    ax.text(0.85, 0.97, "Regulators (RBPs / miRNAs)", ha="center", va="center", fontsize=10, fontweight="bold", color="#475569")

    ax.set_xlim(-0.02, 1.15)
    ax.set_ylim(0.02, 1.02)
    ax.set_title(title, fontsize=12, fontweight="bold", pad=14)


def render_celltype_contrast_network(
    scores_df: pd.DataFrame,
    pop_a: str,
    pop_b: str,
    output_path: Path,
    label_a: str | None = None,
    label_b: str | None = None,
    top_n_edges: int = 30,
    pop_col: str = "population",
    si_mapping: dict[str, float] | None = None,
    dataset_title: str = "scPolASeq Regulatory Rewiring",
):
    """Render a comprehensive 3-panel comparative network: Pop A, Pop B, and Differential Graph."""
    lbl_a = label_a or pop_a
    lbl_b = label_b or pop_b

    edges_a, edges_b, edges_diff = compute_contrast_edge_set(
        scores_df, pop_a, pop_b, top_n_edges=top_n_edges, pop_col=pop_col
    )

    fig, axes = plt.subplots(1, 3, figsize=(22, 11))

    # Panel 1: Pop A Network
    _draw_bipartite_network_panel(
        axes[0], edges_a, weight_col=pop_a,
        title=f"A. {lbl_a} Network",
        is_contrast=False, si_mapping=si_mapping
    )

    # Panel 2: Pop B Network
    _draw_bipartite_network_panel(
        axes[1], edges_b, weight_col=pop_b,
        title=f"B. {lbl_b} Network",
        is_contrast=False, si_mapping=si_mapping
    )

    # Panel 3: Differential Rewiring Network (Delta G = G_B - G_A)
    _draw_bipartite_network_panel(
        axes[2], edges_diff, weight_col="delta_W",
        title=f"C. Differential Rewiring (ΔG = {lbl_b} − {lbl_a})",
        is_contrast=True, si_mapping=si_mapping
    )

    # Legends at bottom
    handles_class = [
        patches.Patch(facecolor=CLASS_PALETTE[k], label=k, edgecolor="k", linewidth=0.5)
        for k in ["ARE-BP", "KH_domain", "RRM", "Zinc_knuckle", "hnRNP", "CSD", "CPA_machinery", "miRNA"]
        if k in CLASS_PALETTE
    ]
    handles_diff = [
        patches.Patch(facecolor="#d95f02", label=f"Gained in {lbl_b} (+ΔW)", edgecolor="k", linewidth=0.5),
        patches.Patch(facecolor="#2b83ba", label=f"Lost in {lbl_b} / Specific to {lbl_a} (−ΔW)", edgecolor="k", linewidth=0.5),
    ]

    fig.legend(handles=handles_class, loc="lower center", bbox_to_anchor=(0.35, 0.02),
               ncol=4, title="Regulator Functional Class (Panels A & B)", fontsize=9, title_fontsize=10, frameon=True)
    fig.legend(handles=handles_diff, loc="lower center", bbox_to_anchor=(0.78, 0.02),
               ncol=2, title="Rewiring Dynamics (Panel C)", fontsize=9, title_fontsize=10, frameon=True)

    plt.suptitle(f"{dataset_title} | Cell-Type Specific Regulatory Networks & Differential Rewiring",
                 fontsize=16, fontweight="bold", y=0.98)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()
    logger.info("Saved cell-type contrast network to %s", output_path)


def export_interactive_celltype_html(
    scores_df: pd.DataFrame,
    pop_a: str,
    pop_b: str,
    output_html: Path,
    label_a: str | None = None,
    label_b: str | None = None,
    top_n_edges: int = 35,
    pop_col: str = "population",
    title: str = "Interactive Cell-Type Regulatory Network",
):
    """Export an interactive Plotly HTML network explorer with dropdown toggles for Pop A, Pop B, and Difference."""
    lbl_a = label_a or pop_a
    lbl_b = label_b or pop_b

    edges_a, edges_b, edges_diff = compute_contrast_edge_set(
        scores_df, pop_a, pop_b, top_n_edges=top_n_edges, pop_col=pop_col
    )

    # Combine all unique genes and regulators
    all_genes = sorted(set(edges_a["gene_name"]).union(set(edges_b["gene_name"])).union(set(edges_diff["gene_name"])))
    all_regs = sorted(set(edges_a["regulator"]).union(set(edges_b["regulator"])).union(set(edges_diff["regulator"])))

    pos_genes = {g: (0.1, y) for g, y in zip(all_genes, np.linspace(0.9, 0.1, len(all_genes)))}
    pos_regs = {r: (0.9, y) for r, y in zip(all_regs, np.linspace(0.9, 0.1, len(all_regs)))}

    fig = go.Figure()

    # Helper to build trace data for a given edge set
    def make_traces(edges: pd.DataFrame, weight_col: str, is_diff: bool = False, visible: bool = False):
        edge_x, edge_y, edge_text = [], [], []
        max_w = edges[weight_col].abs().max() if not edges.empty and edges[weight_col].abs().max() > 0 else 1.0

        for _, r in edges.iterrows():
            g = r["gene_name"]
            reg = r["regulator"]
            w = r[weight_col]
            x0, y0 = pos_genes[g]
            x1, y1 = pos_regs[reg]

            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_text.append(f"{g} ↔ {reg} ({r['regulator_class']})<br>Score: {w:.3f}")

        # Lines trace
        line_color = "#d95f02" if is_diff else "#6366f1"
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y, mode="lines",
            line=dict(width=1.8, color=line_color),
            hoverinfo="none", visible=visible, showlegend=False
        )

        # Genes trace
        gene_trace = go.Scatter(
            x=[pos_genes[g][0] for g in all_genes],
            y=[pos_genes[g][1] for g in all_genes],
            mode="markers+text",
            marker=dict(size=18, color="#3b82f6", symbol="square"),
            text=all_genes, textposition="middle left",
            hoverinfo="text", hovertext=[f"Target Gene: {g}" for g in all_genes],
            visible=visible, showlegend=False
        )

        # Regulators trace
        reg_trace = go.Scatter(
            x=[pos_regs[r][0] for r in all_regs],
            y=[pos_regs[r][1] for r in all_regs],
            mode="markers+text",
            marker=dict(size=14, color="#a855f7", symbol="circle"),
            text=all_regs, textposition="middle right",
            hoverinfo="text", hovertext=[f"Regulator: {r}" for r in all_regs],
            visible=visible, showlegend=False
        )

        return [edge_trace, gene_trace, reg_trace]

    traces_a = make_traces(edges_a, weight_col=pop_a, is_diff=False, visible=True)
    traces_b = make_traces(edges_b, weight_col=pop_b, is_diff=False, visible=False)
    traces_diff = make_traces(edges_diff, weight_col="delta_W", is_diff=True, visible=False)

    for t in traces_a + traces_b + traces_diff:
        fig.add_trace(t)

    # 3 options in dropdown
    fig.update_layout(
        title=f"{title}: {lbl_a} vs {lbl_b}",
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.1, y=1.12,
                buttons=[
                    dict(label=f"1. {lbl_a}", method="update",
                         args=[{"visible": [True, True, True, False, False, False, False, False, False]},
                               {"title": f"Network: {lbl_a}"}]),
                    dict(label=f"2. {lbl_b}", method="update",
                         args=[{"visible": [False, False, False, True, True, True, False, False, False]},
                               {"title": f"Network: {lbl_b}"}]),
                    dict(label=f"3. Difference (ΔG)", method="update",
                         args=[{"visible": [False, False, False, False, False, False, True, True, True]},
                               {"title": f"Differential Rewiring (ΔG = {lbl_b} − {lbl_a})"}]),
                ]
            )
        ],
        height=800,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        margin=dict(l=100, r=100, t=100, b=50),
        plot_bgcolor="#ffffff"
    )

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html))
    logger.info("Saved interactive network HTML to %s", output_html)
