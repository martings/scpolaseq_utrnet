# `utrnet-graph/v0.1`

## Purpose

This is a framework-independent output contract for biological regulatory graph evidence. It wraps existing provider exports for Gemini UTRNNet and other consumers. It does not replace `utrnet-input/v0.1`, define a neural-network representation, or require PyTorch.

## Directory and version metadata

`utrnet export-contract --directory <gene-output>` writes a deterministic bundle under `<gene-output>/utrnet-graph-v0.1/`. The bundle contains `graph_contract.json`, `apa_metric_registry.json`, copied CSV tables, `provenance.json` when present, and copied per-population GraphML/Cytoscape files under `graphs/`.

The manifest records the contract version, provider and repository, source commit when available, source input provenance and hashes, sample/gene context, reference file paths where available, generation configuration, evidence sources, typed relationship categories, canonical entity IDs, file sizes and SHA-256 hashes. Unknown species or genome release stays null.

Canonical IDs use `<type>:<percent-encoded-source-id>`; for example `gene:ENSG00000130066`, `population:Mono`, `pas:site%3A1`, `utr:utr%3Aref%3A1`, `region:utr%3Aref%3A1%3AARE%3A2`, and `regulator:ELAVL1`. Source IDs remain unchanged in CSV columns. IDs are sorted in the manifest so entity indexing does not depend on row order.

## Required tables

| Table | Required fields | Identity / role |
| --- | --- | --- |
| `apa_profiles.csv` | `gene_id`, `celltype_corrected`, `site_id`, `umi_sum`, `total_gene_umi`, `pas_proportion`, `n_cells` | One row per gene, population and PAS. Usage comes from pooled Sierra UMIs; zero coverage remains missing. |
| `utr_isoforms.csv` | `isoform_id`, `gene_id`, `site_id`, `sequence_length`, `ambiguity_status`, `model_weight` | Candidate transcript/3′UTR models; unresolved models are retained. |
| `regulatory_sites.csv` | `regulatory_site_id`, `isoform_id`, `site_id`, `regulator`, `regulator_class`, `genomic_start`, `genomic_end`, `accessibility`, `evidence_tier` | Sequence/imported candidate sites with provenance and unmeasured accessibility preserved. |
| `regulator_abundance.csv` | `population`, `regulator`, `mean_umi_per_cell` | Measured RBP mRNA where available; miRNA abundance can remain missing. |
| `interaction_scores.csv` | `population`, `isoform_id`, `site_id`, `regulatory_site_id`, `regulator`, `regulator_class`, `model_weight`, `accessibility`, `evidence_factor`, `evidence_tier`, `abundance_factor`, `interaction_score`, `score_status` | Candidate scores, never experimental occupancy. |

The validator enforces unique identity keys and cross-table references. Additional columns are allowed and retained.

## Optional tables

The bundle preserves all CSV outputs, including `pas_sites.csv`, `cpa_regulatory_edges.csv`, `cpa_expression.csv`, `cpa_correlations.csv`, `evidence.csv`, `population_metrics.csv`, `rewiring.csv`, `site_retention.csv`, `common_region_accessibility.csv`, and any future additional tables. `pas_sites.csv`, CPA relationships and direct evidence rows are optional for a graph bundle. If `pas_sites.csv` is present, the bundle adds `chr` as an alias of `chrom` for consumer compatibility while preserving both fields.

## Relationship and evidence categories

The manifest names typed relations and their epistemic category:

| Relationship | Category | Provider evidence |
| --- | --- | --- |
| Gene → PAS | `structural` | observed reference association |
| PAS → isoform | `annotation_derived` | transcript model; ambiguity and model weight retained |
| Isoform → regulatory region | `structural` | sequence scan or imported region |
| Region → RBP/miRNA | `inferred_candidate` | motif, seed or imported evidence tier |
| Regulator → gene | `inferred_candidate` | expression- and sequence-weighted candidate score |
| Population → gene expression | `experimental_observation` | curated raw UMI where available |
| Population → PAS/isoform usage | `experimental_observation` | Sierra UMI usage |
| CPA factor → gene/APA | `inferred_candidate` | correlation/literature-derived CPA edge |

GraphML retains existing edge `kind`, `evidence`, score and status attributes. A sequence motif match is not measured binding occupancy. Existing provider evidence tiers such as `motif_heuristic`, `computational_prediction`, `experimental_other_context`, `Sierra_UMI`, and `curated_literature_and_correlation` are preserved verbatim.

## Missingness and ambiguity

- Empty/NaN values remain missing; measured numeric zero remains zero.
- `coverage_status`, `score_status`, accessibility nulls and abundance nulls are not filled in the exported CSVs.
- Contract-bundle `ambiguity_status` maps provider-only names to `unique`, `ambiguous_multiple_models`, or `no_compatible_3utr_model`; original values remain in `source_ambiguity_status`.
- Missing miRNA measurements, no compatible 3′UTR model, zero coverage, and context-incompatible evidence must not be represented as measured zero or experimental support.

## Validation and serialization

`utrnet validate-contract --directory <gene-output-or-bundle>` prints a machine-readable JSON report and exits nonzero on required-table, ID, cross-reference, checksum or GraphML errors. `utrnet export-contract` validates an existing legacy directory first, sorts copied CSV rows with stable keys, uses a deterministic manifest key order, and leaves the legacy outputs unchanged. Re-exporting the same source into a second empty destination produces byte-identical tables and manifest when source inputs are unchanged.

The validator can check provider-side consistency and manifest hashes without importing `UTRNNet`, PyTorch or PyG. Gemini's separate adapter test remains its own integration gate.
