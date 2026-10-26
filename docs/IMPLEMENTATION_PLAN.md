# scPolASeq-UTRNet implementation milestones

This plan follows the read-only PBMC 1k audit. All work stays in a standalone `scpolaseq_utrnet` package and consumes frozen upstream outputs; no alignment or PAS discovery rerun is planned.

| Milestone | Reviewable output | Acceptance check |
| --- | --- | --- |
| 0. Inventory (complete) | `DATA_INVENTORY.md` and JSON, source/contract reconciliation, blockers | Core assets, physical schemas, row counts and checksum coverage documented; absent assets explicit. |
| 1. Contract and ingestion | `utrnet-input/v0.1`, exact-key cell/PAS/count tables, CLI `audit` and `ingest`, synthetic fixture | No barcode suffix stripping; unique join keys; per-cell/per-site Sierra grouping conservation; 1,139 annotations and 1,113 curated cells reconciled. |
| 2. Isoform model | Strand-aware terminal exon/PAS resolver, ambiguity flags, sequence extraction, candidate CLI | Plus/minus, overlapping genes, alternative terminal exons and suspicious spans tested; IFITM3 proximal PAS failed annotation review, so SAT1 and CORO1B form the exploratory panel; ZSWIM6 remains flagged until resolved. |
| 3. Evidence and folding | Versioned RBP/ARE and miRNA site resources, explicit prior-only miRNA mode, cached RNAplfold backend | Evidence tiers and cell-context transfer visible; sequence/parameter hash cache; boundary and common-region accessibility tests. |
| 4. Graph and uncertainty | Structural and inferred edge tables, decomposed scores, population comparisons, null/sensitivity framework | Raw usage, site presence, accessibility, evidence and regulator abundance stay separate; no score described as measured occupancy. |
| 5. PBMC demonstration | SAT1 and CORO1B explorers, static figures, Cytoscape JSON, GraphML, provenance | Runs from frozen PBMC files and from a synthetic fixture; exploratory cross-lineage and transcript-model limits shown; absent miRNA measurements labeled. |
| 6. Quality and release | CI, API/CLI docs, seven modular `.agents/skills`, reprocessing backlog and prioritized GitHub issue drafts | All tests and audit gates pass; Reactome FDR count corrected; no large data, credentials or upstream modifications committed. |

## Immediate implementation choices

- Use the cluster Sierra view as the canonical raw observation source for PBMC 1k. Keep group context and cross-check the cell-type view at the exact `(sample, cell, site)` key before accepting it as a redundant view.
- Store raw `site_id`, `gene_id`, `chrom`, barcode and count columns alongside normalized identifiers. Record coordinate conversion and source-row lineage. Flag `NA` genes and duplicated catalog IDs rather than selecting a record silently.
- Treat `adata_celltyped.h5ad` as the curated population source, with explicit exclusion and QC reasons for the 26 annotation cells absent from it. Use its `counts` layer for measured mRNA abundance; keep miRNA abundance missing unless separately supplied.
- Audit IFITM3 before proposing a 3′UTR regulatory switch; its highlighted proximal PAS did not pass the terminal-3′UTR check in the frozen annotation. Use SAT1 and CORO1B as the first exploratory panel, retaining their ambiguity flags. Keep ZSWIM6 in the audit backlog.

## Issue backlog prepared for review

1. **P0 — Exact-key ingestion and duplicate-site policy.** Address catalog duplicate IDs, 10x suffix preservation, sample scoping, `NA` genes, and Sierra conservation.
2. **P0 — Coordinate and terminal-exon audit.** Pin down source coordinate conventions; validate IFITM3 and audit ZSWIM6 span.
3. **P0 — Correct Reactome FDR reporting.** Count `FDR < 0.05`, fix checkpoint and figure labels, and review cached gene mapping.
4. **P1 — Version regulatory evidence.** Choose RBP/ARE and miRNA resources with accession/version/license/provenance fields.
5. **P1 — Implement reproducible accessibility backend.** RNAplfold configuration, sequence-hash cache and common-region comparison.
6. **P1 — Graph scores and sensitivity.** Explain components, uncertainty, null models and no-reprocessing parameter sweep.
7. **P2 — Explorer, exports, CI and agent skills.** Complete user-facing demonstration and package validation.

These are prepared issue titles and scopes; no GitHub issues have been created.
