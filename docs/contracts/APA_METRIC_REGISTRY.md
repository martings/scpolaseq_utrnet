# APA metric compatibility registry v0.1

## Purpose and migration gate

The registry makes source field, definition, orientation, normalization, range, contrast convention, implementation version, availability, and compatibility status explicit. It is compatibility metadata, not a second implementation of downstream APA metrics. The JSON is included as `apa_metric_registry.json` in graph bundles.

Allowed statuses are `definition_verified`, `definition_ambiguous`, and `definition_pending_migration`. Only an explicitly verified entry with known orientation can pass `interpretation_allowed`. Unknown metric names do not resolve by substring. The registry never renames `si` to DRI and never alters imported values.

## Initial registry

| Metric ID | Source field | Definition / orientation | Status |
| --- | --- | --- | --- |
| `pas_usage_proportion` | `pas_proportion` | `umi_sum / total_gene_umi`; higher means more usage of that PAS; [0,1] | verified |
| `multisite_shortening_index` | `si` | Canonical `sum_i p_i * (i-1)/(N-1)`; proximal = 0, distal = 1; higher = more distal usage | verified against UTRNNet registry and audited downstream code |
| `proximal_pas_fraction` | derived `p_proximal` | Fraction assigned to the first PAS in strand-aware proximal-to-distal order; 1 = proximal | definition verified; no standalone downstream output field |
| `rank_weighted_usage_variant` | not a distinct provider field | Formula and normalization unresolved | ambiguous |
| `fixed_basis_genomic_magnetization` | `magnetization_genomic`, `m_genomic_fixed_a/b`, `delta_m_genomic_fixed` | Strand-aware coordinate weighting on a declared basis; [−1,1], proximal-low/distal-high | verified for audited implementation; basis scope recorded |
| `apa_entropy` | `entropy` | Shannon `−sum p_i log2(p_i)` over PAS passing the configured support threshold; bits | verified for audited implementation |
| `normalized_apa_entropy` | `entropy_norm` | Shannon entropy / `log2(max(2, n_sites_pass))` | verified for audited implementation |
| `differential_apa_statistic` | `delta_si`, `delta_m_genomic_fixed`, `adj_p_value` | Audited `apa.py` stores group A − group B; older tables remain source-specific | pending artifact-level migration |
| `pooled_umi_usage_interval` | `usage_ci_low`, `usage_ci_high` | Jeffreys/Dirichlet marginal interval over pooled PAS UMIs; not biological-replicate uncertainty | verified |

## Strand, order and contrast handling

`order_pas_by_proximity` only orders explicit PAS coordinates on one explicit strand: ascending coordinate on `+`, descending on `-`. The audited downstream `compute_si_and_magnetization` uses the canonical SI formula over PAS meeting the configured per-site UMI threshold (defaults: 5 UMIs/site, 2 passing sites/gene, 10 total UMIs/gene). Thus the support set is part of provenance. The separate `p_proximal` metric is derivable from the first strand-aware PAS and its declared usage denominator, but the audited downstream tables do not emit a standalone column for it. Tests cover two PAS and three-PAS ordering on each strand.

`contrast_delta` returns `group_b - group_a` for verified values when the caller explicitly supplies that order. The audited downstream fields `delta_si` and `delta_m_genomic_fixed` instead store group A − group B; both labels are present in the table. Their existing values are preserved. For those rows, a negative `delta_si` means A has a lower canonical SI (more proximal usage) than B. Historical tables without matching implementation/provenance remain gated and are never sign-inverted.

The fixed-basis genomic metric maps strand-aware PAS coordinates to [−1,1]. The per-group `magnetization_genomic` uses a gene-wide reference span across passing sites, while pairwise `m_genomic_fixed_a/b` uses the union of passing PAS from that pair. Those columns are related but not interchangeable when the bases differ. The rank magnetization in the audited implementation is `2*SI−1`. Shannon entropy and its normalization are verified from the code, but neither is a shortening direction metric.

## Migration path

The local audit compared the registry against `scPolASeq_downstream@365c0722758f7d3a08e9a06050c974207374cfe0` and the UTRNNet metric declarations provided by its maintainer. This registry records those semantics without importing the downstream schema or changing any historical values. Older output tables still require artifact-level implementation and group-order provenance before directional interpretation. The authoritative UTRNNet schema version can be ingested when it is published. Do not overwrite existing PBMC exploratory reports during migration.
