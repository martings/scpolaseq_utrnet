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
        "orientation": "higher_is_more_distal_usage_under_referenced_definition",
        "normalization_basis": "PAS rank from proximal (0) to distal (N-1), not genomic coordinate alone",
        "range": [0.0, 1.0],
        "contrast_convention": "legacy delta fields are retained verbatim; group order is not assumed",
        "source_implementation_version": "not established for imported legacy values",
        "data_availability": "optional imported apa_profiles.si or downstream metric tables",
        "compatibility_status": "definition_pending_migration",
        "reference": "Fabrizio Corzini thesis; downstream source implementation and historical direction audit pending",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "proximal_pas_usage",
        "source_field": "not_available_as_a_distinct_provider_field",
        "definition": "usage assigned to the biologically proximal PAS after strand-aware ordering",
        "orientation": "higher_means_more_proximal_usage",
        "normalization_basis": "total gene PAS UMI within the same population",
        "range": [0.0, 1.0],
        "contrast_convention": "group_b minus group_a only when both groups share a verified PAS ordering",
        "source_implementation_version": "not implemented in provider",
        "data_availability": "not_available",
        "compatibility_status": "definition_pending_migration",
        "reference": "design placeholder; requires authoritative downstream definition",
        "interpretation_allowed": False,
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
        "source_field": "magnetization_genomic",
        "definition": "legacy fixed-basis genomic magnetization; exact basis and historical implementation are not verified here",
        "orientation": "unknown",
        "normalization_basis": "pending authoritative definition",
        "range": None,
        "contrast_convention": "delta_m_genomic_fixed is preserved verbatim; sign is not interpreted",
        "source_implementation_version": "not established for imported legacy values",
        "data_availability": "optional imported apa_profiles.csv or downstream metric tables",
        "compatibility_status": "definition_pending_migration",
        "reference": "downstream metric audit pending",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "apa_entropy",
        "source_field": "entropy",
        "definition": "legacy APA entropy; log base, support, and normalization are not verified",
        "orientation": "unknown",
        "normalization_basis": "pending authoritative definition",
        "range": None,
        "contrast_convention": "no contrast interpretation until source definition is verified",
        "source_implementation_version": "not established for imported legacy values",
        "data_availability": "optional imported apa_profiles.csv or downstream metric tables",
        "compatibility_status": "definition_ambiguous",
        "reference": "downstream metric audit pending",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "normalized_apa_entropy",
        "source_field": "entropy_normalized",
        "definition": "legacy normalized APA entropy; denominator and PAS support are not verified",
        "orientation": "unknown",
        "normalization_basis": "pending authoritative definition",
        "range": None,
        "contrast_convention": "no contrast interpretation until source definition is verified",
        "source_implementation_version": "not established for imported legacy values",
        "data_availability": "optional downstream metric table",
        "compatibility_status": "definition_pending_migration",
        "reference": "downstream metric audit pending",
        "interpretation_allowed": False,
    },
    {
        "metric_id": "differential_apa_statistic",
        "source_field": "delta_si|delta_m_genomic_fixed|adj_p_value",
        "definition": "imported differential APA estimates and support fields; contrast model and group order depend on the source table",
        "orientation": "source_defined; no default direction",
        "normalization_basis": "source-defined",
        "range": None,
        "contrast_convention": "preserve source group_a/group_b and original delta field; never infer direction from name",
        "source_implementation_version": "not established for imported legacy values",
        "data_availability": "optional apa_priority table",
        "compatibility_status": "definition_pending_migration",
        "reference": "downstream metric audit pending",
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


def contrast_delta(metric_id: str, group_a: float, group_b: float) -> float:
    """Compute an explicitly registered group_b-minus-group_a contrast only."""
    entry = metric_definition(metric_id)
    if not interpretation_allowed(metric_id):
        raise ValueError(f"metric {metric_id!r} is not cleared for new contrast interpretation")
    if "group_b minus group_a" not in entry["contrast_convention"]:
        raise ValueError(f"metric {metric_id!r} has no registered contrast convention")
    return float(group_b) - float(group_a)


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
