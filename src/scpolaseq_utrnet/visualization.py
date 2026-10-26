"""Interactive HTML and static publication-style overview."""
from __future__ import annotations

from pathlib import Path
import json
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def _display_regulators(scores: pd.DataFrame, sites: pd.DataFrame, limit: int = 8) -> list[str]:
    """Choose a legible set of regulators without repeating identical site patterns."""
    if scores.empty:
        return []
    ranked = scores.groupby("regulator").interaction_score.sum().sort_values(ascending=False)
    kinds = sites.groupby("regulator").regulator_class.first()
    signatures = {}
    for regulator, group in sites.groupby("regulator"):
        positions = tuple(sorted(set(zip(group.chrom.astype(str), group.genomic_start.astype(int),
                                         group.genomic_end.astype(int), group.strand.astype(str)))))
        signatures[regulator] = (kinds[regulator], positions)
    rbps = [name for name in ranked.index if kinds.get(name) != "miRNA"]
    order = rbps[:1] + [name for name in ranked.index if name not in rbps[:1]]
    selected, seen = [], set()
    for name in order:
        signature = signatures[name]
        if signature in seen:
            continue
        selected.append(name)
        seen.add(signature)
        if len(selected) == limit:
            break
    return selected


def render_gene(directory: Path) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    profiles = pd.read_csv(directory / "apa_profiles.csv")
    scores = pd.read_csv(directory / "interaction_scores.csv")
    isoforms = pd.read_csv(directory / "utr_isoforms.csv")
    pas = pd.read_csv(directory / "pas_sites.csv")
    sites = pd.read_csv(directory / "regulatory_sites.csv")
    retention = pd.read_csv(directory / "site_retention.csv", keep_default_na=False)
    displayed = _display_regulators(scores, sites)
    gene = str(pas.gene_name.iloc[0])
    ordered = pas.sort_values("start", ascending=str(pas.strand.iloc[0]) == "-")
    site_order = ordered.site_id.tolist()
    matrix = profiles.pivot(index="celltype_corrected", columns="site_id", values="pas_proportion")
    matrix = matrix.reindex(columns=site_order).fillna(0)
    fig = make_subplots(rows=2, cols=1, subplot_titles=("Observed PAS usage by population", "Candidate regulatory score"),
                        vertical_spacing=0.18)
    fig.add_trace(go.Heatmap(z=matrix.values, x=[f"PAS {i+1}" for i in range(len(site_order))],
                             y=matrix.index.tolist(), colorscale="Blues", colorbar=dict(title="PAS fraction"),
                             customdata=[[site_order[j] for j in range(len(site_order))] for _ in range(len(matrix))],
                             hovertemplate="%{y}<br>%{customdata}<br>usage=%{z:.3f}<extra></extra>"), row=1, col=1)
    if not scores.empty:
        summary = scores[scores.regulator.isin(displayed)].groupby(["population", "regulator"], as_index=False).interaction_score.sum()
        for regulator, group in summary.groupby("regulator"):
            fig.add_trace(go.Bar(x=group.population, y=group.interaction_score,
                                 name=regulator, hovertemplate="%{x}<br>partial candidate score=%{y:.4f}<extra></extra>"), row=2, col=1)
    fig.update_layout(height=950, title=f"{gene}: observed APA and candidate UTR interactions",
                      barmode="group", margin=dict(l=70, r=40, t=90, b=180))
    fig.update_xaxes(tickangle=-45, row=2, col=1)
    graph_files = sorted((directory / "graphs").glob("*.graphml"))
    network_fig = go.Figure()
    graph_names = []
    kind_colors = {"population": "#111827", "gene": "#2563eb", "PAS": "#f59e0b",
                   "UTR_isoform": "#10b981", "regulatory_region": "#e11d48",
                   "ARE-BP": "#7c3aed", "miRNA": "#9333ea"}
    for file in graph_files:
        graph = nx.read_graphml(file)
        label = graph.graph.get("population", file.stem)
        graph_names.append(label)
        positions = nx.spring_layout(graph, seed=7)
        edge_x, edge_y = [], []
        for source, target in graph.edges():
            a, b = positions[source], positions[target]
            edge_x.extend([a[0], b[0], None])
            edge_y.extend([a[1], b[1], None])
        visible = len(graph_names) == 1
        network_fig.add_trace(go.Scatter(x=edge_x, y=edge_y, mode="lines",
            line=dict(color="#cbd5e1", width=1), hoverinfo="skip", showlegend=False, visible=visible))
        node_ids = list(graph.nodes())
        network_fig.add_trace(go.Scatter(x=[positions[n][0] for n in node_ids],
            y=[positions[n][1] for n in node_ids], mode="markers",
            marker=dict(size=[12 if graph.nodes[n].get("kind") in {"gene", "population"} else 7 for n in node_ids],
                        color=[kind_colors.get(graph.nodes[n].get("kind"), "#64748b") for n in node_ids]),
            text=[f"{graph.nodes[n].get('kind')}: {graph.nodes[n].get('label', n)}<br>"
                  f"{graph.nodes[n].get('ambiguity', '')}" for n in node_ids],
            hovertemplate="%{text}<extra></extra>", showlegend=False, visible=visible))
    if graph_names:
        buttons = []
        for index, label in enumerate(graph_names):
            visible = [False] * (2 * len(graph_names))
            visible[index * 2:index * 2 + 2] = [True, True]
            buttons.append(dict(label=label, method="update", args=[{"visible": visible},
                           {"title": f"Candidate UTR graph: {label}"}]))
        network_fig.update_layout(title=f"Candidate UTR graph: {graph_names[0]}",
            updatemenus=[dict(buttons=buttons, x=0, y=1.15)], height=700,
            xaxis=dict(visible=False), yaxis=dict(visible=False),
            margin=dict(l=20, r=20, t=100, b=20))
    summary = isoforms.ambiguity_status.value_counts().to_dict()
    html = directory / "explorer.html"
    note = ("<h1>scPolASeq-UTRNet: " + gene + "</h1><p><strong>Exploratory model.</strong> "
            "PAS usage is observed in curated PBMC cells. Motifs, accessibility and interactions are predictions "
            "or missing measurements; scores are not molecular occupancy. The sample has no independent "
            "biological replication. miRNA abundance is not inferred from 10x mRNA data.</p>"
            f"<p>Isoform status: {json.dumps(summary)}. Hover over the heatmap to see original PAS IDs. "
            "The bar chart shows up to eight representative regulators with distinct predicted site patterns; "
            "the table and network retain all candidates.</p>")
    site_table = scores[["population", "regulator", "genomic_region", "relative_isoform_usage",
                         "evidence_tier", "accessibility_status", "abundance_source",
                         "interaction_score", "score_ci_low", "score_ci_high", "score_status",
                         "population_coverage_status"]].to_html(index=False, escape=True) if not scores.empty else "<p>No candidate motif interactions for this gene.</p>"
    html.write_text("<!doctype html><html><head><meta charset='utf-8'><title>UTRNet explorer</title>"
                    "<style>body{font-family:system-ui,sans-serif;max-width:1500px;margin:auto;padding:24px}"
                    "table{border-collapse:collapse;font-size:12px}th,td{padding:5px;border:1px solid #ddd}</style>"
                    "</head><body>" + note + fig.to_html(full_html=False, include_plotlyjs=True) +
                    network_fig.to_html(full_html=False, include_plotlyjs=False) +
                    "<h2>Site retention across candidate PAS models</h2><div style='overflow:auto;max-height:350px'>" +
                    retention.to_html(index=False, escape=True) + "</div>" +
                    "<h2>Decomposed candidate interactions</h2><div style='overflow:auto;max-height:500px'>" +
                    site_table + "</div></body></html>")

    static = directory / "summary.png"
    fig2, axes = plt.subplots(2, 1, figsize=(15, max(8, len(matrix) * 0.35 + 4)), constrained_layout=True)
    axes[0].imshow(matrix.values, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    axes[0].set_yticks(range(len(matrix.index)), matrix.index, fontsize=7)
    axes[0].set_xticks(range(len(site_order)), [str(i + 1) for i in range(len(site_order))], fontsize=7)
    axes[0].set_title("Observed PAS proportions (5′→3′ in transcript order)")
    if not scores.empty:
        plot = scores[scores.regulator.isin(displayed)].groupby(["population", "regulator"]).interaction_score.sum().unstack(fill_value=0)
        plot.plot(kind="bar", ax=axes[1], legend=False)
        axes[1].legend(title="Representative candidates", bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        axes[1].set_ylabel("Candidate score (partial if inputs missing)")
        axes[1].set_title("Up to eight distinct predicted site patterns; full set in explorer")
    else:
        axes[1].text(0.5, 0.5, "No sequence motif interactions passed annotation checks",
                     transform=axes[1].transAxes, ha="center")
    axes[1].tick_params(axis="x", labelrotation=45)
    fig2.savefig(static, dpi=200)
    plt.close(fig2)
    return html, static
