---
name: utr-graph-builder
description: Build population-specific heterogeneous UTR graphs with separate structural and candidate interaction edges and decomposed scores.
---

# UTR graph builder

Input: observed APA profiles, candidate isoforms, regulatory sites, optional accessibility, regulator abundance and evidence tiers.

Keep gene, PAS, isoform, motif/region, regulator and population nodes distinct. Structural edges describe annotation or presence; interaction edges are candidate predictions. Export every score factor and missingness status. Do not call scores molecular occupancy. Compare populations by genomic regulatory region and regulator, with gained/lost availability and sensitivity to sparse populations.

Output: score tables, population graphs, rewiring tables, Cytoscape JSON and GraphML. Check graph joins, per-PAS usage normalization, ambiguity weights, export round-trip and cell coverage. Handoff scientific claims to QA rather than labeling them validated.
