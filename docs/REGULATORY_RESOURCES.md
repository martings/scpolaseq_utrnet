# Regulatory resource provenance

## RBP motifs

The two bundled motif records are a minimal extract of the [official ATtRACT download](https://attract.cnic.es/download), whose [documentation](https://attract.cnic.es/documentation) describes experimentally supported RBP-associated motifs and their PubMed and experiment fields. Archive URL: `https://attract.cnic.es/attract/static/ATtRACT.zip`. Archive SHA-256: `adf73651e8e4d75c0fea86a5fb10926c93b7b04628d065789c1196b9ee0a3a76`. The archive was last modified 2016-09-27 according to its HTTP headers; `2016-09-27` is the local resource version label. `scripts/prepare_attract.py` verifies the archive checksum and regenerates `resources/are_motifs.tsv` exactly.

| RBP | Motif | ATtRACT matrix | PubMed | Source assay/context |
| --- | --- | --- | --- | --- |
| ELAVL1 | AUUUA | 1103 | 9882309 | UV cross-linking and EMSA; HeLa nuclear extracts/recombinant protein. |
| ZFP36 | AUUUA | 38 | 18256032 | Immunoprecipitation-based evidence; no PBMC-specific occupancy. |

The experimentally supported relationship is RBP ↔ motif in the source study. Each SAT1 motif location is a **new computational sequence match**. Its graph edge is therefore a candidate interaction transferred from another context. The `evidence.csv` output retains assay, PubMed ID, matrix ID and source URL.

## miRNA sites

The PBMC configuration scans the bundled [MirGeneDB 3.0 human mature miRNA FASTA](https://www.mirgenedb.org/fasta/hsa?mat=1). The 23,323-byte download has SHA-256 `d1499ef63043414d57165c26fb95858c9547eca1ded54aacbbdc1385bda6f322` and contains 559 mature sequences. For each sequence, the scanner takes miRNA positions 2–8 and finds exact reverse-complement 7mer-m8 matches in reconstructed, spliced 3′UTRs. It emits isoform and genomic coordinates only when all seven genomic bases are contiguous. These are **computational sequence priors**, not experimentally validated miRNA–target relationships, and no PBMC miRNA expression or AGO occupancy was measured. Seed-only matching has many false positives; graph scores remain partial where miRNA abundance or RNAplfold accessibility is missing.

The separate optional genomic-site TSV accepts `gene_id`, `regulator`, `chrom`, 1-based inclusive `start/end`, `strand`, `evidence_tier`, `source`, `version`, and `url`. This can hold site-resolved external evidence after assembly and coordinate validation. The [miRDB download page](https://www.mirdb.org/download.html) offers predicted gene targets; gene-level predictions alone cannot establish which reconstructed PAS isoform retains a specific site. The synthetic fixture demonstrates the TSV import and marks abundance as missing-prior-only.
