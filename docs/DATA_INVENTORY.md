# PBMC 1k data inventory — Phase 0

Audit date: 2026-10-25 UTC. This is a read-only audit of the earlier frozen `pbmc_1k_v3` run. The companion [machine-readable inventory](DATA_INVENTORY.json) lists 210 physical files with portable `${UPSTREAM}`/`${DOWNSTREAM}` paths, sizes, row counts where applicable, SHA-256 checksums for 152 files, explicit checksum omissions for large/binary files, and preliminary completeness flags. The current MVP uses the completed 2026-10-24 replica, recorded separately in [PBMC1K_REPLICA_INVENTORY.json](PBMC1K_REPLICA_INVENTORY.json) and [PBMC1K_VALIDATION.md](PBMC1K_VALIDATION.md). Source tables have no declared schema version; they are recorded as `unversioned`. No upstream computation was run.

## Source locations and validated assets

| Asset | Physical source | Finding |
| --- | --- | --- |
| PAS site catalog | `${UPSTREAM}/apa/site_catalog.tsv` | 821,840 rows; 38,068 repeated `site_id` rows, with no conflicting gene IDs among duplicates. Includes interval coordinates and mixed `chr` prefixes. |
| PAS reference | `${UPSTREAM}/pas_reference/pas_reference.tsv` | 741,392 unique `site_id` values; every ID occurs in the catalog, while 42,380 catalog IDs are absent here. `start == end` in every reference row, so it is a point representation. Build manifest reports 25 bp merge distance. |
| Sierra counts | `${UPSTREAM}/sierra_quant/{cluster,cell_type}/` | Cluster and cell-type views each contain 3,054,180 rows and 7,779,859 UMIs. All rows have a known annotation barcode and PAS reference ID; no zero or negative counts. The cell-type view is `unlabeled`. Treat views as alternative groupings of the same observations. |
| Cell identifiers | `labels/cell_annotations.tsv`, `labels/group_map.tsv`, `provenance/pbmc_1k_v3.barcode_registry.tsv` | 1,139 rows in each, 1,139 distinct canonical cell IDs and corrected barcodes; registry and annotation barcodes match. `group_map` contains only `cluster`, not curated cell types. |
| Expression | `labels/pbmc_1k_v3.clustered.h5ad`; downstream `qc/adata_qc_processed.h5ad` and `cell_labels/adata_celltyped.h5ad` | Shapes: 1,137 × 2,000; 1,137 × 2,000; 1,113 × 61,228. Only curated downstream AnnData has a `counts` layer and `celltype_corrected`. Its 1,113 cell IDs are all in the 1,139-row annotations; 26 annotated cells are not in the curated object. |
| Downstream APA | `outputs/pbmc1k/apa/` | 482,666 site-support rows; 22,130 SI/magnetization/entropy profiles; 52,750 comparison/priority rows. A second 52,750-row fixed-basis table exists under `apafield/`. These are derived summaries, not independent observations. |
| Reference | `reference/reference.genome.fa`, `.fai`, `reference.annotation.gtf`, `terminal_exons.tsv` | FASTA, FAI and GTF are present; GTF header identifies GRCh38, GENCODE v44 / Ensembl 110. Terminal exon table has 252,835 rows. The reference manifest points to work-directory paths rather than the published file paths; provenance should use both recorded and resolved paths. |
| Source code | Downstream package `src/scpolaseq_downstream/` | `io.py`, `apa.py`, `celltyping.py`, `validation.py`, `enrichment.py`, and plotting/config modules are available. Reuse formulas and audit logic where appropriate, but preserve original IDs in the new ingestion contract. |

## Compatibility with upstream documentation

Compared with upstream [outputs.md](https://github.com/martings/scPolASeq/blob/main/docs/outputs.md) and [module_contracts.md](https://github.com/martings/scPolASeq/blob/main/docs/module_contracts.md):

- `reference/reference_manifest.tsv` is documented but the run has `reference/reference_manifest.json`.
- The documented alignment-local barcode registry and alignment manifest are published under `provenance/` in this run.
- The documented grouped coverage count/summary TSVs, `apa/apa_usage.tsv`, `apa/site_usage_matrix.tsv`, `apa/gene_summary.tsv`, `apa/apa_stats.tsv`, and `pas_scoring/pas_scored_events.tsv` are absent. The run instead has `apa_calls/` event and PDUI tables. These absent files are not assumed to be available.
- The label contract lists `cell_type`, but its physical values are all `unlabeled`. Curated types are downstream in `adata_celltyped.h5ad`; do not infer them from `group_map` or the Sierra `cell_type/unlabeled` export.
- The sample sheet has two library FASTQ rows, `L001` and `L002`, merged into one sample-scoped analysis ID. Physical annotation and Sierra `library_id` are `pbmc_1k_v3`; preserve original sample-sheet library IDs separately rather than inventing per-cell library attribution.
- Site catalog `chrom` uses both `X`/`11` and `chrX`/`chr11`. PAS reference uses `chr` names. Catalog intervals and reference point coordinates must retain their original convention and be normalized explicitly against the FASTA index.
- The 18 `*.bw` coverage files are zero bytes. BEDGraph tracks exist; BigWig availability is false.

## Validation findings and blockers

1. **Join cardinality and unmapped genes.** The site catalog contains 235,892 rows with literal gene ID `NA`; the PAS reference contains 234,526. The downstream priority table has 24 `NA` gene records, including spans of 247,127,995 bp. Exclude unresolved gene IDs from candidate graphs until mapped or flagged with provenance. A many-to-one site catalog join requires deduplication with an explicit rule; a raw merge expands 741,392 PAS records to 760,196 rows.
2. **Coordinate model.** PAS reference points are inside corresponding catalog intervals after matching site IDs, but the source formats differ. A matching `site_id` is not proof of an unambiguous transcript or terminal exon assignment. The BED files are zero-based half-open; TSV point and interval conventions require confirmation in source code before reconstruction.
3. **Barcode semantics.** The upstream docs describe possible 10x `-1` suffixes, but this physical run has none in its 1,139 annotation barcodes or sampled STARsolo barcodes. Downstream `normalize_barcode()` strips prefixes and suffixes and joins by normalized barcode alone. The new layer must join by sample, library, and exact corrected barcode, and fail on collisions.
4. **Count conservation.** Cluster and cell-type Sierra views match in total rows and UMIs, but this aggregate check does not establish per-cell/per-site equality. Implement exact key-level conservation in the ingestion milestone before publication.
5. **Candidate annotation.** IFITM3 (`ENSG00000142089`) is a plausible minus-strand technical candidate with multiple sites and terminal exons, pending transcript/PAS compatibility checks. ZSWIM6 (`ENSG00000130449`) has a 210,297 bp fixed-basis span in the priority table; some catalog sites lie well upstream of the one terminal exon observed in this table. Hold it out of the initial validated panel pending a genomic mapping audit.
6. **Statistical interpretation.** The Reactome output has 20 displayed pathways, minimum FDR 0.923 and zero pathways with FDR < 0.05. The notebook checkpoint records `n_reactome_pathways_fdr05: 20` by counting displayed rows, which is incorrect. The 942 `robust` APA comparison rows are cross-lineage exploratory results from one biological sampling context, not independent biological replication.
7. **Regulator evidence.** No versioned RBP motif library, miRNA target database, small-RNA abundance, or RNAplfold output was found in the audited output trees. Mark these as missing until an evidence resource is chosen and versioned; never substitute poly(A)-selected 10x mRNA counts for endogenous miRNA abundance.

## Reusable downstream logic and cautions

`apa.py` computes SI from observed site proportions and ordered PAS rank; magnetization and entropy formulas can be adopted with tests. `io.py` streams Sierra TSVs and derives gene IDs, but its barcode normalization is lossy and its site-to-cell mapping depends on the current one-sample context. The new contract should import only safe helpers or adapt formulas with explicit key validation. `enrichment.py` and the notebook Reactome summary need a corrected FDR filter before any pathway claim.

## Decision for the next milestone

The frozen output set is sufficient to start a versioned, exact-ID ingestion layer and a limited IFITM3 annotation check. It is not yet sufficient to claim a validated 3′UTR regulatory network: transcript/PAS assignment, regulatory evidence resources, miRNA priors, and accessibility predictions remain to be established. See [IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) for reviewable milestones and acceptance checks.
