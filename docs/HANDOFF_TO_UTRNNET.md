# Handoff: `utrnet-graph/v0.1` to Gemini UTRNNet

## Contract version and artifacts

Provider version: **`utrnet-graph/v0.1`**. The provider writes a portable bundle at `<gene-output>/utrnet-graph-v0.1/`. It contains the contract manifest, metric registry, required/available CSVs, optional provenance and copied GraphML/Cytoscape artifacts. Current root-level CSVs remain available for older consumers.

The provider-side command interface is:

```bash
utrnet build-graph --config configs/pbmc1k.yaml --gene SAT1
utrnet validate-contract --directory outputs/pbmc1k/ENSG00000130066/utrnet-graph-v0.1
```

For a previously generated directory:

```bash
utrnet export-contract --directory outputs/pbmc1k/ENSG00000130066
utrnet validate-contract --directory outputs/pbmc1k/ENSG00000130066/utrnet-graph-v0.1
```

No library import is required by a consumer to inspect the CSV bundle.

## Required tables and exact consumer fields

| File | Required fields |
| --- | --- |
| `apa_profiles.csv` | `gene_id`, `celltype_corrected`, `site_id`, `umi_sum`, `total_gene_umi`, `pas_proportion`, `n_cells` |
| `utr_isoforms.csv` | `isoform_id`, `gene_id`, `site_id`, `sequence_length`, `ambiguity_status`, `model_weight` |
| `regulatory_sites.csv` | `regulatory_site_id`, `isoform_id`, `site_id`, `regulator`, `regulator_class`, `genomic_start`, `genomic_end`, `accessibility`, `evidence_tier` |
| `regulator_abundance.csv` | `population`, `regulator`, `mean_umi_per_cell` |
| `interaction_scores.csv` | `population`, `isoform_id`, `site_id`, `regulatory_site_id`, `regulator`, `regulator_class`, `model_weight`, `accessibility`, `evidence_factor`, `evidence_tier`, `abundance_factor`, `interaction_score`, `score_status` |

All other columns are retained. Optional `pas_sites.csv`, CPA tables and evidence/provenance tables are preserved when available. A present PAS table receives a `chr` alias for its current `chrom` field without removing the source field.

## Guarantees and mapping notes

- The provider contract uses stable, reversible IDs and sorts the contract copy deterministically. Original source identifiers remain in table fields.
- The manifest contains provider commit when discoverable, generation configuration, input-file references/hashes, evidence sources, file hashes, entity IDs, table/GraphML row counts and relationship evidence categories.
- Structural, annotation-derived, experimental-observation and inferred-candidate evidence remain distinguishable. Candidate scores are not binding occupancy.
- Empty values remain missing; measured zeros remain zero. Coverage, accessibility, abundance, unresolved models and score statuses survive the copy.
- Provider ambiguity values `ambiguous_sequence_models` and `same_sequence_multiple_transcripts` map to `ambiguous_multiple_models` in the bundle, with originals retained in `source_ambiguity_status`. `unique` and `no_compatible_3utr_model` remain explicit.
- Relationship categories distinguish gene→PAS, PAS→isoform, isoform→region, region→regulator, regulator→gene, population→expression/usage, and optional CPA edges.
- `apa_metric_registry.json` ships with the bundle. `si`, magnetization, entropy and legacy deltas are preserved, but pending/ambiguous records are not eligible for new directional biological claims.

## How Gemini can consume it

UTRNNet's current `cell-regulatory-state/v0.1` can build its `HeteroData` from the five required tables. Keep the existing graph conversion and anti-leakage policy in UTRNNet. For the new provider contract, point the reader at the bundle directory and use its stable entity inventory and metric registry in addition to the core CSVs.

The inspected Gemini feature branch is `feat/m0-m2-utrnnet-hetero-gnn`. It already contains `src/utrnnet/contracts/utrnet_graph_v0_1.py`, `src/utrnnet/graphs/adapter.py`, `src/utrnnet/contracts/schema.py`, and `tests/gnn/test_upstream_contract_compatibility.py`. Provider-side tests in this PR cover deterministic export and validation; the Gemini owner should extend the compatibility test to validate the bundle path and verify entity/relationship counts plus missingness and ambiguity round-trip.

Provider tests can be run independently:

```bash
python -m pytest -q tests/test_graph_contract.py tests/test_metric_registry.py tests/test_signature.py
```

The provider branch does not execute UTRNNet or modify its active feature branch. The current Gemini adapter imports helpers from `scpolaseq_utrnet.cis`; its owner should remove that private package import before UTRNNet claims full independence from provider internals.

## Unresolved metric semantics and limitations

The cited thesis formula associates higher multisite SI with more distal usage when PAS order is proximal-to-distal. This provider has not migrated historical dataset definitions or contrast group order, so `si` stays `definition_pending_migration` and cannot drive validated shortening/lengthening claims. Genomic magnetization, entropy variants and differential APA fields also remain pending or ambiguous. The active downstream metric audit is the authority for a future migration.

The provider bundle is a biological graph evidence export, not a Tahoe matrix and not a PyTorch dataset. It does not add timepoints, gene perturbations, measured miRNA abundance, unsupported pathway activity, or causal target claims. Non-present modalities remain explicit in `regulatory-signature/v0.1`.
