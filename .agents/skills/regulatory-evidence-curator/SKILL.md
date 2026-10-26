---
name: regulatory-evidence-curator
description: Version and annotate RBP, ARE and miRNA regulatory-site resources for UTRNet without conflating motifs with binding evidence.
---

# Regulatory evidence curator

Input: reconstructed isoforms and candidate motif/target databases with version, accession, license, assay and context metadata.

Record each site in isoform and genomic coordinates, regulator class, evidence tier, source URL and whether evidence is PBMC-specific or transferred. Treat generic ARE patterns as motif hypotheses. Distinguish imported experimental binding from computational prediction. Standard 10x poly(A)-selected mRNA data do not measure endogenous miRNA abundance; use labeled prior-only mode until small-RNA data exist.

Output: regulatory-site and evidence tables with reproducible resource versions. Check coordinate retention, duplicate evidence, motif null/background comparison and missing provenance. Handoff only evidence whose origin and limits are inspectable.
