# APA metric compatibility registry v0.1

## Purpose and migration gate

The registry makes source field, definition, orientation, normalization, range, contrast convention, implementation version, availability, and compatibility status explicit. It is compatibility metadata, not a second implementation of downstream APA metrics. The JSON is included as `apa_metric_registry.json` in graph bundles.

Allowed statuses are `definition_verified`, `definition_ambiguous`, and `definition_pending_migration`. Only an explicitly verified entry with known orientation can pass `interpretation_allowed`. Unknown metric names do not resolve by substring. The registry never renames `si` to DRI and never alters imported values.

## Initial registry

| Metric ID | Source field | Definition / orientation | Status |
| --- | --- | --- | --- |
| `pas_usage_proportion` | `pas_proportion` | `umi_sum / total_gene_umi`; higher means more usage of that PAS; [0,1] | verified |
| `multisite_shortening_index` | `si` | Reference formula: `sum_i p_i * (i-1)/(N-1)` with PAS proximal-to-distal; referenced orientation is higher = more distal usage | pending migration |
| `proximal_pas_usage` | not a distinct provider field | Requires verified strand-aware proximal PAS selection; higher = more proximal usage | pending migration |
| `rank_weighted_usage_variant` | not a distinct provider field | Formula and normalization unresolved | ambiguous |
| `fixed_basis_genomic_magnetization` | `magnetization_genomic` | Legacy fixed-basis value preserved; orientation and basis unresolved | pending migration |
| `apa_entropy` | `entropy` | Log base, support and normalization unresolved | ambiguous |
| `normalized_apa_entropy` | `entropy_normalized` | Denominator and PAS support unresolved | pending migration |
| `differential_apa_statistic` | `delta_si`, `delta_m_genomic_fixed`, `adj_p_value` | Preserve source group order and statistic; no default direction | pending migration |
| `pooled_umi_usage_interval` | `usage_ci_low`, `usage_ci_high` | Jeffreys/Dirichlet marginal interval over pooled PAS UMIs; not biological-replicate uncertainty | verified |

## Strand, order and contrast handling

`order_pas_by_proximity` only orders explicit PAS coordinates on one explicit strand: ascending coordinate on `+`, descending on `-`. It does not calculate SI, change historical ranks, or make imported metric values authoritative. Tests cover two PAS and three-PAS ordering on each strand.

`contrast_delta` returns `group_b - group_a` only for registry entries explicitly verified with that contrast convention. It rejects legacy SI and magnetization until group ordering and historical direction are reconciled. A negative contrast therefore means only that the supplied B value is numerically lower; biological shortening/lengthening is not inferred for pending metrics.

## Migration path

When the downstream audit publishes an authoritative definition, update the relevant registry entries with its source field, exact implementation version, range and contrast order. Mark the entry `definition_verified` only after fixture comparisons against the audited implementation pass. Keep original fields and provenance in backward-compatible exports. Do not overwrite existing PBMC exploratory reports during migration.
