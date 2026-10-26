"""Stage invalidation planner; default dry-run performs no computation."""
from __future__ import annotations

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
