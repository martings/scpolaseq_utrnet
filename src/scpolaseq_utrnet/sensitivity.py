"""Stage invalidation planner and graph-only frozen-score sensitivity sweep."""
from __future__ import annotations

import numpy as np
import pandas as pd

PARAMETERS = {
    "score_weights": ("graph_only", "graph scores, rewiring, explorer", "seconds", False),
    "display_threshold": ("graph_only", "explorer", "seconds", False),
    "motif_library": ("sequence", "regulatory sites, scores, graphs", "minutes", False),
    "rna_window": ("sequence", "fold cache, scores, graphs", "minutes to hours", False),
    "min_population_cells": ("downstream_APA", "profiles, scores, graphs", "minutes", False),
    "min_pas_umi": ("downstream_APA", "profiles, scores, graphs", "minutes", False),
    "pas_merge_distance": ("upstream_PAS", "PAS reference, Sierra, all downstream outputs", "hours", True),
    "scapture_activation": ("upstream_PAS", "PAS catalog, PAS reference, Sierra, all downstream outputs", "hours to days", True),
}


def plan(changes: list[str]) -> list[dict]:
    result = []
    for name in changes:
        if name not in PARAMETERS:
            raise ValueError(f"unknown sensitivity parameter {name}")
        group, invalidates, runtime, upstream = PARAMETERS[name]
        result.append({"parameter": name, "group": group, "invalidates": invalidates,
                       "expected_runtime": runtime, "upstream_reprocessing": upstream,
                       "validation": "compare gene/PAS membership, coverage, rankings and graph changes"})
    return result


def graph_sweep(scores: pd.DataFrame, min_population_cells: int = 15) -> pd.DataFrame:
    """Re-score frozen candidate edges without modifying PAS or sequence inputs.

    This is a model-assumption check, not a biological replication analysis.
    Missing abundance remains a neutral *partial* factor in every scenario.
    """
    columns = ["scenario", "population", "included", "n_cells", "candidate_edges",
               "score_sum", "score_delta_from_baseline", "rank_correlation",
               "top_ten_overlap", "miRNA_evidence_multiplier", "missing_accessibility_factor",
               "min_population_cells"]
    if scores.empty:
        return pd.DataFrame(columns=columns)
    required = {"population", "regulator", "genomic_region", "regulator_class", "relative_isoform_usage", "site_present",
                "model_weight", "accessibility", "evidence_factor", "abundance_factor",
                "interaction_score", "n_cells"}
    if not required.issubset(scores.columns):
        raise ValueError(f"scores lack {sorted(required - set(scores.columns))}")
    if min_population_cells < 1:
        raise ValueError("min_population_cells must be positive")
    scores = scores.reset_index(drop=True)
    scenarios = [
        ("baseline", 1.0, 1.0, min_population_cells),
        ("miRNA_evidence_half", 0.5, 1.0, min_population_cells),
        ("miRNA_evidence_double", 2.0, 1.0, min_population_cells),
        ("missing_accessibility_half", 1.0, 0.5, min_population_cells),
        ("stricter_cell_coverage", 1.0, 1.0, max(25, min_population_cells + 10)),
    ]
    usage = scores.relative_isoform_usage.fillna(0).to_numpy(dtype=float)
    common = (usage * scores.site_present.to_numpy(dtype=float) *
              scores.model_weight.to_numpy(dtype=float) *
              scores.evidence_factor.to_numpy(dtype=float) *
              scores.abundance_factor.fillna(1).to_numpy(dtype=float))
    baseline = scores.interaction_score.to_numpy(dtype=float)
    rebuilt = common * scores.accessibility.fillna(1).to_numpy(dtype=float)
    if not np.allclose(rebuilt, baseline, rtol=1e-8, atol=1e-10, equal_nan=True):
        raise ValueError("frozen score factors do not reconstruct the baseline")
    mirna = scores.regulator_class.eq("miRNA").to_numpy()
    rows = []
    for name, mirna_multiplier, missing_accessibility, threshold in scenarios:
        accessibility = scores.accessibility.fillna(missing_accessibility).to_numpy(dtype=float)
        values = common * accessibility * np.where(mirna, mirna_multiplier, 1.0)
        if mirna_multiplier == 1.0 and missing_accessibility == 1.0:
            values = baseline.copy()
        for population, group in scores.groupby("population", sort=True):
            positions = group.index.to_numpy()
            n_cells = int(group.n_cells.iloc[0])
            included = n_cells >= threshold
            base = baseline[positions]
            candidate = values[positions]
            edges = group[["regulator", "genomic_region"]].copy()
            edges["baseline"] = base
            edges["candidate"] = candidate
            edges = edges.groupby(["regulator", "genomic_region"], sort=True).sum()
            correlation = (edges.baseline.rank().corr(edges.candidate.rank())
                           if included and len(edges) > 1 else np.nan)
            top_count = min(10, len(edges))
            top_base = set(edges.sort_values("baseline", ascending=False, kind="stable").head(top_count).index)
            top_candidate = set(edges.sort_values("candidate", ascending=False, kind="stable").head(top_count).index)
            rows.append({"scenario": name, "population": population, "included": included,
                         "n_cells": n_cells, "candidate_edges": len(edges) if included else 0,
                         "score_sum": float(candidate.sum()) if included else np.nan,
                         "score_delta_from_baseline": float(candidate.sum() - base.sum()) if included else np.nan,
                         "rank_correlation": correlation,
                         "top_ten_overlap": len(top_base & top_candidate) / top_count if included else np.nan,
                         "miRNA_evidence_multiplier": mirna_multiplier,
                         "missing_accessibility_factor": missing_accessibility,
                         "min_population_cells": threshold})
    return pd.DataFrame(rows, columns=columns)
