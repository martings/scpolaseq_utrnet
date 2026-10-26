"""Heterogeneous graph construction and decomposed interaction scores."""
from __future__ import annotations

import json
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd


EVIDENCE_FACTOR = {"motif_heuristic": 0.25, "computational_prediction": 0.35,
                   "external_database": 0.55, "experimental_other_context": 0.75,
                   "experimental_matching_context": 1.0}


def site_retention(isoforms: pd.DataFrame, regulatory_sites: pd.DataFrame) -> pd.DataFrame:
    """Summarize which candidate PAS models retain each predicted genomic site."""
    columns = ["regulator", "regulator_class", "genomic_region", "present_pas_ids",
               "absent_pas_ids", "model_dependent_pas_ids", "unresolved_pas_ids",
               "retention_status"]
    if regulatory_sites.empty:
        return pd.DataFrame(columns=columns)
    compatible = isoforms.loc[isoforms.sequence_length > 0]
    models = compatible.groupby("site_id").isoform_id.apply(set).to_dict()
    unresolved = sorted(set(isoforms.site_id) - set(models))
    frame = regulatory_sites.copy()
    frame["genomic_region"] = (frame.chrom.astype(str) + ":" + frame.genomic_start.astype(str) +
                               "-" + frame.genomic_end.astype(str) + ":" + frame.strand.astype(str))
    rows = []
    for (regulator, kind, region), group in frame.groupby(
            ["regulator", "regulator_class", "genomic_region"], sort=True):
        present, absent, dependent = [], [], []
        observed = set(group.isoform_id)
        for pas_id, candidates in sorted(models.items()):
            overlap = candidates & observed
            if not overlap:
                absent.append(pas_id)
            elif overlap == candidates:
                present.append(pas_id)
            else:
                dependent.append(pas_id)
        status = ("model_dependent" if dependent else
                  "pas_specific" if absent else "retained_across_compatible_pas")
        rows.append({"regulator": regulator, "regulator_class": kind, "genomic_region": region,
                     "present_pas_ids": ";".join(present), "absent_pas_ids": ";".join(absent),
                     "model_dependent_pas_ids": ";".join(dependent),
                     "unresolved_pas_ids": ";".join(unresolved), "retention_status": status})
    return pd.DataFrame(rows, columns=columns)


def population_metrics(scores: pd.DataFrame) -> pd.DataFrame:
    """Descriptive candidate-network size and weight concentration by population."""
    columns = ["population", "n_cells", "population_coverage_status", "candidate_rows", "candidate_edges",
               "candidate_regulators", "candidate_genomic_regions", "total_candidate_score",
               "score_concentration_hhi", "partial_or_missing_rows"]
    if scores.empty:
        return pd.DataFrame(columns=columns)
    rows = []
    for population, group in scores.groupby("population", sort=True):
        weights = group.groupby(["regulator", "genomic_region"]).interaction_score.sum().clip(lower=0)
        total = float(weights.sum())
        rows.append({"population": population, "n_cells": int(group.n_cells.iloc[0]),
                     "population_coverage_status": group.population_coverage_status.iloc[0],
                     "candidate_rows": len(group), "candidate_edges": len(weights),
                     "candidate_regulators": group.regulator.nunique(),
                     "candidate_genomic_regions": group.genomic_region.nunique(),
                     "total_candidate_score": total,
                     "score_concentration_hhi": float(((weights / total) ** 2).sum()) if total else float("nan"),
                     "partial_or_missing_rows": int((group.score_status != "complete_prediction").sum())})
    return pd.DataFrame(rows, columns=columns)


def score_interactions(profiles: pd.DataFrame, isoforms: pd.DataFrame,
                       regulatory_sites: pd.DataFrame, abundance: pd.DataFrame,
                       min_population_cells: int) -> pd.DataFrame:
    columns = ["population", "gene_id", "isoform_id", "site_id", "regulatory_site_id",
               "regulator", "regulator_class", "genomic_region", "relative_isoform_usage",
               "site_present", "model_weight", "accessibility", "accessibility_status",
               "evidence_factor", "evidence_tier", "abundance_factor", "abundance_source",
               "interaction_score", "score_ci_low", "score_ci_high", "score_status",
               "population_coverage_status", "n_cells", "gene_umi"]
    if regulatory_sites.empty:
        return pd.DataFrame(columns=columns)
    usage = profiles.rename(columns={"celltype_corrected": "population", "pas_proportion": "relative_isoform_usage"})
    frame = regulatory_sites.merge(isoforms[["isoform_id", "model_weight", "gene_id"]], on="isoform_id", validate="many_to_one")
    frame = frame.merge(usage[["population", "site_id", "relative_isoform_usage", "usage_ci_low", "usage_ci_high", "n_cells", "total_gene_umi"]],
                        on="site_id", validate="many_to_many")
    frame = frame.copy()
    frame["gene_umi"] = frame.total_gene_umi
    frame["site_present"] = 1.0
    frame["evidence_factor"] = frame.evidence_tier.map(EVIDENCE_FACTOR)
    if frame.evidence_factor.isna().any():
        raise ValueError("unrecognized evidence tier")
    frame = frame.merge(abundance[["regulator", "population", "mean_umi_per_cell", "abundance_source"]],
                        on=["regulator", "population"], validate="many_to_one")
    maxima = frame.groupby("regulator").mean_umi_per_cell.transform("max")
    frame["abundance_factor"] = np.where(frame.mean_umi_per_cell.notna(),
                                          (frame.mean_umi_per_cell + 0.1) / (maxima + 0.1), np.nan)
    frame["genomic_region"] = (frame.chrom.astype(str) + ":" + frame.genomic_start.astype(str) +
                                "-" + frame.genomic_end.astype(str) + ":" + frame.strand.astype(str))
    # Neutral factors preserve a partial ranking; score_status prevents treating it as complete evidence.
    factor = (frame.site_present * frame.model_weight * frame.accessibility.fillna(1) *
              frame.evidence_factor * frame.abundance_factor.fillna(1))
    frame["interaction_score"] = frame.relative_isoform_usage.fillna(0) * factor
    frame["score_ci_low"] = frame.usage_ci_low * factor
    frame["score_ci_high"] = frame.usage_ci_high * factor
    frame["score_status"] = np.where(frame.accessibility.notna() & frame.abundance_factor.notna(),
        "complete_prediction", "partial_missing_accessibility_or_abundance")
    frame.loc[frame.relative_isoform_usage.isna(), "score_status"] = "no_observed_gene_umi"
    frame["population_coverage_status"] = np.where(frame.n_cells >= min_population_cells,
                                                     "meets_min_cells", "small_population_exploratory")
    return frame[columns]


def build_graph(scores: pd.DataFrame, profiles: pd.DataFrame, isoforms: pd.DataFrame,
                regulatory_sites: pd.DataFrame, population: str,
                cpa_edges: pd.DataFrame | None = None,
                cpa_expr: pd.DataFrame | None = None) -> nx.MultiDiGraph:
    graph = nx.MultiDiGraph(population=population, model="utrnet-v0.2")
    p = profiles[profiles.celltype_corrected == population]
    s = scores[scores.population == population]
    if p.empty:
        raise ValueError(f"unknown population {population}")
    gene = str(p.gene_id.iloc[0])
    pop_node = f"population:{population}"
    gene_node = f"gene:{gene}"
    pop_attrs = {"kind": "population", "label": population}
    if "si" in p.columns and pd.notna(p["si"].iloc[0]):
        pop_attrs["shortening_index"] = float(p["si"].iloc[0])
    if "magnetization_genomic" in p.columns and pd.notna(p["magnetization_genomic"].iloc[0]):
        pop_attrs["magnetization_genomic"] = float(p["magnetization_genomic"].iloc[0])
    if "entropy" in p.columns and pd.notna(p["entropy"].iloc[0]):
        pop_attrs["entropy"] = float(p["entropy"].iloc[0])
    graph.add_node(pop_node, **pop_attrs)
    graph.add_node(gene_node, kind="gene", label=gene)
    for row in p.itertuples(index=False):
        pas_id = f"pas:{row.site_id}"
        graph.add_node(pas_id, kind="PAS", label=row.site_id)
        graph.add_edge(gene_node, pas_id, kind="gene_has_PAS", evidence="observed_reference")
        for iso in isoforms[isoforms.site_id == row.site_id].itertuples(index=False):
            graph.add_node(iso.isoform_id, kind="UTR_isoform", label=iso.isoform_id,
                           ambiguity=iso.ambiguity_status, span_status=iso.span_status,
                           sequence_length=int(iso.sequence_length))
            graph.add_edge(pas_id, iso.isoform_id, kind="PAS_defines_candidate_isoform", evidence="annotation_model")
            graph.add_edge(pop_node, iso.isoform_id, kind="uses_isoform", evidence="Sierra_UMI",
                           weight=float(row.pas_proportion) if pd.notna(row.pas_proportion) else 0.0,
                           umi=int(row.umi_sum))
    for row in s.itertuples(index=False):
        motif_node = f"region:{row.regulatory_site_id}"
        regulator_node = f"regulator:{row.regulator}"
        graph.add_node(motif_node, kind="regulatory_region", label=row.genomic_region)
        graph.add_node(regulator_node, kind=row.regulator_class, label=row.regulator)
        graph.add_edge(row.isoform_id, motif_node, kind="contains_site", evidence="sequence_or_import")
        graph.add_edge(motif_node, regulator_node, kind="candidate_binding", evidence=row.evidence_tier,
                       weight=float(row.interaction_score), score_status=row.score_status)
        graph.add_edge(regulator_node, gene_node, kind="candidate_regulation", evidence=row.evidence_tier,
                       weight=float(row.interaction_score), score_status=row.score_status)

    if cpa_edges is not None and not cpa_edges.empty:
        pop_cpa_expr = {}
        if cpa_expr is not None and not cpa_expr.empty:
            sub_expr = cpa_expr[cpa_expr.population == population]
            for r in sub_expr.itertuples(index=False):
                pop_cpa_expr[r.factor] = r

        for row in cpa_edges.itertuples(index=False):
            cpa_node = str(row.source)
            target_node = gene_node if row.target == "gene:target" else str(row.target)
            if target_node in graph or target_node == gene_node:
                node_attrs = {
                    "kind": "cpa_factor",
                    "label": row.factor,
                    "complex": row.complex,
                    "pmid": str(row.pmid),
                }
                if row.factor in pop_cpa_expr:
                    r_exp = pop_cpa_expr[row.factor]
                    node_attrs["mean_expr"] = float(r_exp.mean_expr)
                    node_attrs["pct_expressing"] = float(r_exp.pct_expressing)

                graph.add_node(cpa_node, **node_attrs)
                graph.add_edge(
                    cpa_node,
                    target_node,
                    kind=row.kind,
                    weight=float(row.correlation_r),
                    p_value=float(row.p_value),
                    direction=row.direction,
                    evidence=row.evidence_tier,
                )
    return graph


def compare_populations(scores: pd.DataFrame, priority_df: pd.DataFrame | None = None) -> pd.DataFrame:
    if scores.empty:
        return pd.DataFrame(columns=["population_a", "population_b", "regulator", "genomic_region",
                                     "score_a", "score_b", "delta_score", "availability_change"])
    summary = scores.groupby(["population", "regulator", "genomic_region"], as_index=False).agg(
        interaction_score=("interaction_score", "sum"),
        score_ci_low=("score_ci_low", "sum"), score_ci_high=("score_ci_high", "sum"),
        available=("relative_isoform_usage", lambda x: bool((x > 0).any())))
    populations = sorted(summary.population.unique())
    records = []

    # Map priority comparisons if available
    prio_lookup = {}
    if priority_df is not None and not priority_df.empty:
        for r in priority_df.itertuples(index=False):
            key1 = (getattr(r, "group_a", None), getattr(r, "group_b", None))
            key2 = (getattr(r, "group_b", None), getattr(r, "group_a", None))
            prio_lookup[key1] = r
            prio_lookup[key2] = r

    for i, a in enumerate(populations):
        for b in populations[i + 1:]:
            left = summary[summary.population == a].set_index(["regulator", "genomic_region"])
            right = summary[summary.population == b].set_index(["regulator", "genomic_region"])
            prio_row = prio_lookup.get((a, b))

            for key in left.index.union(right.index):
                la, rb = left.loc[key] if key in left.index else None, right.loc[key] if key in right.index else None
                sa = float(la.interaction_score) if la is not None else 0.0
                sb = float(rb.interaction_score) if rb is not None else 0.0
                available_a = bool(la.available) if la is not None else False
                available_b = bool(rb.available) if rb is not None else False
                rec = {
                    "population_a": a, "population_b": b, "regulator": key[0],
                    "genomic_region": key[1], "score_a": sa, "score_b": sb,
                    "delta_score": sb - sa,
                    "delta_ci_low": (float(rb.score_ci_low) if rb is not None else 0.0) -
                                    (float(la.score_ci_high) if la is not None else 0.0),
                    "delta_ci_high": (float(rb.score_ci_high) if rb is not None else 0.0) -
                                     (float(la.score_ci_low) if la is not None else 0.0),
                    "availability_change": "gained" if available_b and not available_a else
                    "lost" if available_a and not available_b else "retained",
                }
                if prio_row is not None:
                    rec["delta_m_genomic_fixed"] = getattr(prio_row, "delta_m_genomic_fixed", None)
                    rec["delta_si"] = getattr(prio_row, "delta_si", None)
                    rec["endpoint_artifact_suspect"] = getattr(prio_row, "endpoint_artifact_suspect", False)
                    rec["is_robust"] = getattr(prio_row, "is_robust", False)
                    rec["priority_score"] = getattr(prio_row, "priority_score", None)
                records.append(rec)
    return pd.DataFrame(records)


def export_graphs(scores: pd.DataFrame, profiles: pd.DataFrame, isoforms: pd.DataFrame,
                  regulatory_sites: pd.DataFrame, directory: Path,
                  cpa_edges: pd.DataFrame | None = None,
                  cpa_expr: pd.DataFrame | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for population in sorted(profiles.celltype_corrected.unique()):
        graph = build_graph(scores, profiles, isoforms, regulatory_sites, population,
                            cpa_edges=cpa_edges, cpa_expr=cpa_expr)
        slug = "".join(c.lower() if c.isalnum() else "_" for c in population).strip("_")
        (directory / f"{slug}.cytoscape.json").write_text(json.dumps(nx.cytoscape_data(graph), indent=2))
        nx.write_graphml(graph, directory / f"{slug}.graphml")
