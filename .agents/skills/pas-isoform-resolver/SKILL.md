---
name: pas-isoform-resolver
description: Resolve PAS points against strand-aware GTF terminal exons and reconstruct candidate 3′UTR sequences with ambiguity flags.
---

# PAS isoform resolver

Input: versioned PAS reference, GTF, indexed genome FASTA, and gene candidates. Preserve original coordinates and chromosome names alongside any normalized representation.

Require a PAS to overlap a compatible terminal exon and transcript-oriented 3′UTR. Build spliced sequence through the PAS on the correct strand. Keep distinct transcript models, alternative terminal exons and ambiguous sequence associations visible. Flag suspicious spans and sites without compatible models; never force a two-isoform interpretation.

Output: isoform table, sequence coordinate map and per-site status. Check plus/minus orientation, overlap, sequence length, genomic bounds and alternative exon cases. Handoff only candidates whose status and limits are documented; send unresolved spans to scientific QA.
