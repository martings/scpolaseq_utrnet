---
name: scpolaseq-output-auditor
description: Audit frozen scPolASeq run products, physical schemas, barcodes, counts, and provenance before UTRNet ingestion.
---

# scPolASeq output auditor

Input: a run directory, downstream output directory, and source contract documents. Read files and manifests without changing upstream results.

Produce a human and machine-readable inventory with paths, checksums or explicit checksum gaps, row counts, columns, completeness and validation findings. Reconcile documented outputs with actual files. Check exact sample/library/barcode keys, PAS ID joins, coordinate conventions, Sierra grouping redundancy and AnnData layers. Do not infer a field from documentation alone.

Handoff when the frozen input set, unresolved joins and missing resources are explicit. Send ambiguous mappings to the isoform resolver and count discrepancies to scientific QA. Do not start upstream computation as part of an audit.
