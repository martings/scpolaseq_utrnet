"""scPolASeq-UTRNet public API."""

from .config import Config, load_config
from .pipeline import build_gene_graphs
from .cpa import CPA_MACHINERY_CATALOG, extract_cpa_expression, correlate_cpa_with_apa
from .cerna import export_cerna_analysis
from .dynamics import compute_graph_velocity
from .cis import (
    annotate_motif_survival_probabilities,
    compute_fractional_cis_dosage,
    compute_dynamic_dosage_flux,
)
from .flow import (
    compute_edge_flows,
    compute_node_divergence,
    compare_flow_divergence_trajectories,
)
from .graph_viz import (
    filter_top_regulatory_edges,
    compute_contrast_edge_set,
    render_celltype_contrast_network,
    export_interactive_celltype_html,
)

__all__ = [
    "Config",
    "load_config",
    "build_gene_graphs",
    "CPA_MACHINERY_CATALOG",
    "extract_cpa_expression",
    "correlate_cpa_with_apa",
    "export_cerna_analysis",
    "compute_graph_velocity",
    "annotate_motif_survival_probabilities",
    "compute_fractional_cis_dosage",
    "compute_dynamic_dosage_flux",
    "compute_edge_flows",
    "compute_node_divergence",
    "compare_flow_divergence_trajectories",
    "filter_top_regulatory_edges",
    "compute_contrast_edge_set",
    "render_celltype_contrast_network",
    "export_interactive_celltype_html",
]
__version__ = "0.2.2"
