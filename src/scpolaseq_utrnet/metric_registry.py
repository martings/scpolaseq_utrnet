"""Versioned APA metric metadata; legacy numeric values are never rewritten."""
from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterable, Mapping

REGISTRY_VERSION = "apa-metric-registry/v0.1"

# Historical fields are intentionally kept distinct from metric identifiers.
# SI and magnetization semantics remain gated until the downstream audit lands.
METRIC_REGISTRY: list[dict[str, Any]] = [
    {
        "metric_id": "pas_usage_proportion",
        "source_field": "pas_proportion",
        "definition": "umi_sum / total_gene_umi for one PAS in one population; null when total_gene_umi is zero",
        "orientation": "higher_means_greater_usage_of_this_pas",
        "normalization_basis": "pooled curated-cell Sierra UMIs across PAS for the gene and population",
        "range": [0.0, 1.0],
        "contrast_convention": "difference is group_b minus group_a when explicitly requested",
        "source_implementation_version": "scpolaseq_utrnet/0.2.0",
        "data_availability": "apa_profiles.csv",
        "compatibility_status": "definition_verified",
        "reference": "src/scpolaseq_utrnet/ingest.py:make_profiles",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "multisite_shortening_index",
        "source_field": "si",
        "definition": "sum_i p_i * (i - 1) / (N - 1), with PAS ordered proximal to distal",
        "orientation": "distal_high_proximal_low; 0=proximal, 1=distal",
        "normalization_basis": "Within-group UMI proportions over PAS passing min_umi_per_site; ranks run proximal (0) to distal (N-1). Current downstream defaults: min_umi_per_site=5, min_sites_per_gene=2, min_total_umi_per_gene=10; configuration is per run.",
        "range": [0.0, 1.0],
        "contrast_convention": "historical downstream delta_si is group_a minus group_b; new contrasts require an explicit orientation and group labels. Preserve legacy deltas verbatim.",
        "source_implementation_version": "scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0 audited; older output provenance is dataset-specific",
        "data_availability": "si_and_magnetization_by_celltype.tsv::si; celltype_specific_APA_priority_table.tsv::si_a,si_b,delta_si",
        "compatibility_status": "definition_verified",
        "reference": "scPolASeq_downstream/src/scpolaseq_downstream/apa.py::compute_si_and_magnetization and _test_pair_shared_basis; UTRNNet APA Metric Registry",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "proximal_pas_fraction",
        "source_field": "p_proximal (derived; not emitted as a standalone downstream field)",
        "definition": "the PAS usage fraction assigned to the first PAS in strand-aware proximal-to-distal order",
        "orientation": "proximal_high_distal_low; 1=proximal",
        "normalization_basis": "Within-observation gene PAS UMI proportions using the same declared PAS support set and denominator; retain the selected proximal site ID and support threshold in provenance.",
        "range": [0.0, 1.0],
        "contrast_convention": "explicit group_b minus group_a for newly derived contrasts; source legacy contrast fields retain their own recorded order",
        "source_implementation_version": "UTRNNet APA Metric Registry; derived metric, not emitted by audited downstream apa.py",
        "data_availability": "derivable from strand-aware PAS ordering plus site UMI proportions; no standalone p_proximal output in audited downstream tables",
        "compatibility_status": "definition_verified",
        "reference": "UTRNNet APA Metric Registry; scPolASeq_downstream/src/scpolaseq_downstream/apa.py::parse_pas_coordinates",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "rank_weighted_usage_variant",
        "source_field": "not_available_as_a_distinct_provider_field",
        "definition": "rank-weighted PAS usage variant; formula not yet reconciled",
        "orientation": "unknown",
        "normalization_basis": "pending authoritative definition",
        "range": None,
        "contrast_convention": "pending authoritative definition",
        "source_implementation_version": "unknown",
        "data_availability": "not_available",
        "compatibility_status": "definition_ambiguous",
        "reference": "downstream metric audit pending",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "fixed_basis_genomic_magnetization",
        "source_field": "magnetization_genomic|m_genomic_fixed_a|m_genomic_fixed_b|delta_m_genomic_fixed",
        "definition": "sum_i p_i * d_i on a shared strand-aware genomic PAS coordinate basis, where d_i maps the proximal-to-distal coordinate range linearly to [-1, 1]; zero when the basis span is zero",
        "orientation": "distal_high_proximal_low",
        "normalization_basis": "Single-group magnetization uses the gene-wide coordinate span across PAS passing the threshold in any group. Pairwise fields use the union of PAS passing min_umi_per_site in group_a or group_b as the shared basis; p_i is each group's usage over its passing counts, with absent basis PAS assigned zero. Pairwise M is therefore not interchangeable with the single-group field when bases differ.",
        "range": [-1.0, 1.0],
        "contrast_convention": "downstream delta_m_genomic_fixed is m_genomic_fixed_a minus m_genomic_fixed_b; group labels are explicit. Preserve legacy deltas verbatim.",
        "source_implementation_version": "scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0 audited; older output provenance is dataset-specific",
        "data_availability": "si_and_magnetization_by_celltype.tsv::magnetization_genomic; celltype_specific_APA_priority_table.tsv::m_genomic_fixed_a,m_genomic_fixed_b,delta_m_genomic_fixed",
        "compatibility_status": "definition_verified",
        "reference": "scPolASeq_downstream/src/scpolaseq_downstream/apa.py::compute_si_and_magnetization and _test_pair_shared_basis; UTRNNet APA Metric Registry",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "apa_entropy",
        "source_field": "entropy",
        "definition": "Shannon entropy -sum_i p_i * log2(p_i) over PAS passing the configured per-site UMI threshold",
        "orientation": "higher_means_more_distributed_usage; not directional shortening",
        "normalization_basis": "Within-group UMI proportions over passing PAS; entropy is measured in bits",
        "range": None,
        "contrast_convention": "no contrast interpretation until source definition is verified",
        "source_implementation_version": "scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0 audited; older output provenance is dataset-specific",
        "data_availability": "si_and_magnetization_by_celltype.tsv::entropy",
        "compatibility_status": "definition_verified",
        "reference": "scPolASeq_downstream/src/scpolaseq_downstream/apa.py::compute_si_and_magnetization",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "normalized_apa_entropy",
        "source_field": "entropy_norm",
        "definition": "Shannon entropy divided by log2(max(2, n_sites_pass))",
        "orientation": "higher_means_more_distributed_usage; not directional shortening",
        "normalization_basis": "Normalized over the number of PAS passing the configured per-site UMI threshold; intended range [0,1] for normalized proportions",
        "range": [0.0, 1.0],
        "contrast_convention": "explicit group contrast only; source table's group order must be retained",
        "source_implementation_version": "scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0 audited; older output provenance is dataset-specific",
        "data_availability": "si_and_magnetization_by_celltype.tsv::entropy_norm",
        "compatibility_status": "definition_verified",
        "reference": "scPolASeq_downstream/src/scpolaseq_downstream/apa.py::compute_si_and_magnetization",
        "interpretation_allowed": True,
    },
    {
        "metric_id": "differential_apa_statistic",
        "source_field": "delta_si|delta_m_genomic_fixed|adj_p_value",
        "definition": "In the audited apa.py implementation, delta_si=si_a-si_b and delta_m_genomic_fixed=m_genomic_fixed_a-m_genomic_fixed_b; group_a/group_b are stored. Other historical tables may use different implementations and require artifact-level provenance.",
        "orientation": "distal_high_proximal_low for SI and fixed-basis genomic M; delta sign follows group_a minus group_b",
        "normalization_basis": "SI uses each group's passing PAS ranks; fixed-basis M uses the shared pairwise PAS coordinate basis",
        "range": None,
        "contrast_convention": "audited output fields are group_a minus group_b; preserve source group_a/group_b and original deltas. Other legacy tables remain source-defined.",
        "source_implementation_version": "scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0 audited; older output provenance is dataset-specific",
        "data_availability": "celltype_specific_APA_priority_table.tsv",
        "compatibility_status": "definition_pending_migration",
        "reference": "scPolASeq_downstream/src/scpolaseq_downstream/apa.py::_test_pair_shared_basis",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "pooled_umi_usage_interval",
        "source_field": "usage_ci_low|usage_ci_high",
        "definition": "Jeffreys Dirichlet marginal interval from pooled PAS UMI counts",
        "orientation": "interval bounds around pas_proportion",
        "normalization_basis": "pooled gene PAS UMI; not biological-replicate uncertainty",
        "range": [0.0, 1.0],
        "contrast_convention": "interval subtraction is not a replicate-level differential test",
        "source_implementation_version": "scpolaseq_utrnet/0.2.0",
        "data_availability": "apa_profiles.csv when generated",
        "compatibility_status": "definition_verified",
        "reference": "src/scpolaseq_utrnet/ingest.py:make_profiles",
        "interpretation_allowed": True,
    },
]


def metric_definition(metric_id: str) -> dict[str, Any]:
    """Return a copy of a registry record by stable metric ID."""
    for entry in METRIC_REGISTRY:
        if entry["metric_id"] == metric_id:
            return deepcopy(entry)
    raise KeyError(f"unknown APA metric_id: {metric_id}")


def legacy_metric_definition(source_field: str) -> dict[str, Any] | None:
    """Resolve a legacy field exactly; ambiguous or missing names are never guessed."""
    matches = [entry for entry in METRIC_REGISTRY if source_field in entry["source_field"].split("|")]
    if len(matches) != 1:
        return None
    return deepcopy(matches[0])


def interpretation_allowed(metric_id: str) -> bool:
    """Gate new biological direction claims on a verified definition and explicit orientation."""
    entry = metric_definition(metric_id)
    return bool(entry["compatibility_status"] == "definition_verified" and entry["interpretation_allowed"]
                and entry["orientation"] not in {"unknown", "source_defined; no default direction"})


def contrast_delta(metric_id: str, group_a: float, group_b: float, *, orientation: str) -> float:
    """Compute a new contrast only when its group direction is explicitly supplied."""
    entry = metric_definition(metric_id)
    if not interpretation_allowed(metric_id):
        raise ValueError(f"metric {metric_id!r} is not cleared for new contrast interpretation")
    if orientation == "group_b_minus_group_a":
        return float(group_b) - float(group_a)
    if orientation == "group_a_minus_group_b":
        return float(group_a) - float(group_b)
    raise ValueError("orientation must be 'group_a_minus_group_b' or 'group_b_minus_group_a'")


def contrast_record(
    metric_id: str, group_a_label: str, group_a: float, group_b_label: str, group_b: float,
    *, orientation: str, source_implementation: str,
) -> dict[str, Any]:
    """Return a machine-readable contrast with preserved labels and direction provenance."""
    if not group_a_label or not group_b_label or not source_implementation:
        raise ValueError("group labels and source_implementation are required")
    return {
        "metric_id": metric_id,
        "group_a": {"label": group_a_label, "value": float(group_a)},
        "group_b": {"label": group_b_label, "value": float(group_b)},
        "delta": contrast_delta(metric_id, group_a, group_b, orientation=orientation),
        "orientation": orientation,
        "source_implementation": source_implementation,
    }


def order_pas_by_proximity(sites: Iterable[Mapping[str, Any]]) -> list[str]:
    """Order PAS identifiers proximal-to-distal using strand-aware coordinates only.

    This helper defines coordinate ordering; it does not calculate or validate SI.
    """
    rows = list(sites)
    if not rows:
        return []
    strands = {str(row.get("strand", "")) for row in rows}
    if len(strands) != 1 or next(iter(strands)) not in {"+", "-"}:
        raise ValueError("PAS ordering requires one explicit strand per gene")
    strand = next(iter(strands))
    for row in rows:
        if not row.get("site_id") or row.get("start") is None:
            raise ValueError("each PAS requires site_id and start")
    return [str(row["site_id"]) for row in sorted(
        rows, key=lambda row: ((float(row["start"]) if strand == "+" else -float(row["start"])), str(row["site_id"]))
    )]


def registry_document() -> dict[str, Any]:
    return {"registry_version": REGISTRY_VERSION, "metrics": deepcopy(METRIC_REGISTRY)}
