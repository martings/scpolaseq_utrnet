# PBMC 1k MVP validation

## Frozen inputs

Current upstream input: the completed 2026-10-24 PBMC 1k replica, supplied through `UTRNET_UPSTREAM`. Its pinned downstream notebook output, supplied through `UTRNET_DOWNSTREAM`, reports technical concordance: 1,113 QC-pass cells, 16 curated cell types, 22,130 SI profiles, 52,750 differential rows, 942 robust APA comparison rows and 614 unique robust genes. Its Reactome audit correctly reports 0 pathways at FDR ≤ 0.05 among 20 displayed. Scientific review remains pending.

## UTRNet result

On this run, `utrnet build-graph --gene SAT1` produced 26 reference PAS, 4,764 gene-specific Sierra rows, 9 PAS with at least one compatible GTF 3′UTR sequence, 2 distinct AUUUA genomic regions, 16 RBP motif records across transcript models, and 354 exact miRNA 7mer-m8 seed records for 30 distinct mature miRNAs across isoform models. There are 5,920 population-specific candidate interaction rows, 4,560 pairwise rewiring rows, and GraphML/Cytoscape exports for all 16 curated types. The AUUUA associations come from checksum-pinned ATtRACT records for ELAVL1 and ZFP36; their occurrences in SAT1 are sequence predictions. The miRNA seed records come from checksum-pinned MirGeneDB 3.0 mature sequences and are also predictions. All miRNA abundance entries are missing-prior-only; 5,550 interaction rows carry partial-score status and 370 have no observed gene UMI in that population/PAS combination. `utrnet render --gene SAT1` produced an interactive HTML explorer and static PNG. These generated outputs remain in ignored `outputs/` and can be reproduced from the frozen inputs.

The second panel gene, CORO1B, has 10 reference PAS, 6 with compatible GTF 3′UTR sequence, 4 AUUUA motif records, and 876 exact miRNA seed records across isoform models. Its 17 output isoform rows include 12 `ambiguous_sequence_models`, 4 unresolved PAS, and 1 unique model. The 241 distinct regulator–genomic-region site-retention rows include 204 PAS-specific and 37 model-dependent cases. Its 16 population-specific graphs, explorer, and static PNG were generated from the same frozen inputs. The visualization limits its overview bars to eight distinct site patterns while retaining the complete candidate set in tables and graph exports. CORO1B is exploratory because most compatible models are ambiguous; it does not demonstrate experimentally established binding.

The first IFITM3 audit found 11 reference PAS and 374 Sierra rows, but the highlighted proximal PAS at `chr11:320565` has no compatible annotated 3′UTR model; no bundled AUUUA motif edge was generated. This is a reason to validate the annotation before using IFITM3 as a regulatory switch example. ZSWIM6 remains on hold because its reported 210,297 bp PAS span conflicts with a simple terminal-exon interpretation.

## Remaining validation gates

- The completed replica has only a cluster Sierra export. Exact `(sample, cell, PAS)` conservation across grouping views cannot be checked within this run; the audit reports the missing comparison view explicitly. The earlier run had matching aggregate totals, and the synthetic fixture exercises exact-key comparison.
- Review SAT1 transcript/PAS ambiguity and compare genomic motif regions across models before biological interpretation.
- Add PBMC-specific RBP/AGO occupancy or independently validated miRNA sites. ATtRACT motifs and MirGeneDB seed matches do not establish in-cell binding.
- Run RNAplfold with recorded parameters and check common-region accessibility changes; it was not installed for this build.
- Apply uncertainty intervals or null models for rewiring and keep cross-lineage conclusions exploratory.
