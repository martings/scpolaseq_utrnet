# Stable graph contract audit

## Scope and operational boundary

This audit covers the current `scpolaseq_utrnet/main` implementation and the read-only Gemini interface branch `martings/UTRNNet:feat/m0-m2-utrnnet-hetero-gnn`. The provider's PBMC graph builder, CSV tables, GraphML, Cytoscape JSON, and `provenance.json` already work; this change adds a contract bundle and validation around those outputs rather than rebuilding graph construction.

No scPolASeq, Sierra, SCAPTURE, Nextflow, dengue, or GPU workload was run. Synthetic fixtures are the only execution data used by the new contract tests. Work was performed in a temporary isolated clone at `/tmp/codex-utrnet-graph-contract`; the shared `/home/megas/agents/scpolaseq_utrnet` checkout was observed clean and left untouched.

## Existing provider behavior

- `build_gene_graphs` writes gene-scoped `apa_profiles.csv`, `pas_sites.csv`, `utr_isoforms.csv`, `regulatory_sites.csv`, `regulator_abundance.csv`, `interaction_scores.csv`, `evidence.csv`, `provenance.json`, per-population GraphML/Cytoscape JSON, and optional CPA and comparison tables.
- `graphs.build_graph` already encodes gene→PAS, PAS→candidate isoform, population→isoform usage, isoform→regulatory site, region→regulator candidate binding, regulator→gene candidate regulation, and optional CPA relations.
- `interaction_scores.csv` distinguishes `score_status`; `regulatory_sites.csv` carries evidence tier, source, sequence coordinates and possibly missing accessibility. `utr_isoforms.csv` retains unresolved and ambiguous models.
- `utrnet-input/v0.1` is the existing input contract. It is not the output graph contract.

## Consumer compatibility inspected

The Gemini feature branch defines `cell-regulatory-state/v0.1` and an upstream reader for `utrnet-graph/v0.1`. It requires the five tables and fields specified in [the handoff](HANDOFF_TO_UTRNNET.md), derives stable typed node identifiers, preserves missing-value masks, and maps evidence tiers. Its synthetic integration test calls the provider's existing graph builder and checks entity conservation.

Compatibility gaps found and addressed in the provider bundle:

1. Provider ambiguity labels include `ambiguous_sequence_models` and `same_sequence_multiple_transcripts`, while the consumer enumerates `ambiguous_multiple_models`. The contract copy maps both to the canonical category and preserves the original in `source_ambiguity_status`; legacy source CSVs are unchanged.
2. Provider PAS reference exports use `chrom`; the consumer's optional-table schema names `chr`. The bundle keeps `chrom` and adds the alias `chr` when absent.
3. The existing consumer validator checks required columns and basic counts but does not verify table identity uniqueness or cross-table references. The provider validator now checks these at export time.
4. The Gemini adapter branch currently imports `scpolaseq_utrnet.cis` helper functions. Its migration document says this should be decoupled before upstream dependency pinning; this provider PR does not change UTRNNet.

## APA metric compatibility risk

`make_profiles` computes PAS UMI proportions and pooled-UMI Jeffreys/Dirichlet intervals. It also imports legacy `si`, magnetization and entropy fields from downstream outputs. The maintainer-confirmed UTRNNet registry defines canonical SI as `sum_i p_i * (i-1)/(N-1)`, with proximal = 0 and distal = 1. I audited `scPolASeq_downstream/src/scpolaseq_downstream/apa.py` at commit `365c0722758f7d3a08e9a06050c974207374cfe0`: the formula and orientation match.

The audited implementation filters each group's PAS by `min_umi_per_site` before ranking; defaults are 5 UMIs/site, 2 passing sites/gene, and 10 total UMIs/gene. SI then uses group-specific ranks on that passing set. Its pairwise `delta_si` is `si_a - si_b`, with `group_a` and `group_b` written to each row. `delta_m_genomic_fixed` is likewise `m_genomic_fixed_a - m_genomic_fixed_b`. No sign flip or historical table rewrite is part of this PR. Older tables remain gated unless their producing implementation and group order are known.

The fixed-basis genomic metric maps strand-aware PAS order to [−1,1]. Per-group `magnetization_genomic` uses a gene-wide passing-site coordinate span; the pairwise `m_genomic_fixed_a/b` fields use the union of passing sites for that pair, so the values can differ when support sets differ. `magnetization_rank` is `2*SI−1`. `entropy` is Shannon entropy in bits over passing PAS; `entropy_norm` divides by `log2(max(2, n_sites_pass))`. Proximal PAS fraction `p_proximal` is a UTRNNet-defined derived metric, but the audited downstream tables do not emit it as a standalone field. The full definitions and this source-vs-history distinction are in [the APA registry](contracts/APA_METRIC_REGISTRY.md).

## Backward compatibility and remaining limits

- `utrnet-input/v0.1`, current root CSV files, GraphML/Cytoscape files, and existing graph functions remain in place. Contract output is an additive `utrnet-graph-v0.1/` directory.
- The contract validates table schemas, uniqueness, references, GraphML readability, provenance hashes, deterministic entity inventory and copied-file hashes. It does not assert that candidate motifs are experimental binding or that a metric definition is biologically correct.
- Species/reference values remain null in metadata when the source provenance does not provide them; the adapter does not infer them from a gene symbol.
- The provider still has one-gene-at-a-time exports and the PBMC example remains the only real small demonstration. The signature contract is generic but its non-immune fixture is synthetic.
- No external UTRNNet code was imported or executed in this isolated provider test run. Gemini's exact adapter compatibility test must be run in that repository by its owner.
