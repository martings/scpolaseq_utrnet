# `utrnet-input/v0.1` data contract

Input source is a frozen scPolASeq run plus a curated downstream AnnData. Preserve raw identifiers and coordinates; normalization creates additional columns without rewriting originals.

| Table | Primary key / required fields | Source and caveat |
| --- | --- | --- |
| `cells` | `(sample_id, library_id, barcode_corrected)` and reversible `cell_id`; curated cell type, cluster, QC | Upstream annotations left-joined by exact `cell_id` to downstream AnnData. Cells absent from curated AnnData remain flagged. |
| `pas_sites` | `site_id`, stable `pas_reference_id`, gene ID, chrom, 1-based PAS point, strand, reference source | Published PAS reference; site catalog intervals are a separate coordinate representation. |
| `pas_counts` | source file, library ID, group level/ID, raw `cell_barcode`, `site_id`, nonnegative integer UMI, joined `cell_id` | Cluster Sierra view is the canonical observation stream. Cell-type Sierra view represents the same UMIs and must not be summed with it. |
| `gene_expression` | `cell_id`, gene ID, count value, `layer=counts`, `normalization=raw_UMI` | Curated AnnData sparse `counts` layer; mRNA only. |
| `apa_profiles` | gene, population, PAS, UMI, n cells, gene UMI total, PAS proportion, coverage status | Pooled UMI estimator over curated cells; zero-coverage profiles are missing, not zero proportions. |
| `utr_isoforms` | isoform ID, PAS ID, transcript IDs, sequence, terminal exon, model weight, ambiguity status | GTF exon/CDS/UTR and FASTA; distinct sequence alternatives for one PAS share usage by explicit model weight. |
| `regulatory_sites` | stable region ID, isoform/genomic coordinates, regulator, source, evidence tier, accessibility status | Motif scan, mature-miRNA seed scan, or imported miRNA resource. A sequence match is not measured occupancy. |
| `site_retention` | regulator and genomic region; present, absent, model-dependent and unresolved PAS IDs | Site availability across candidate 3′UTR models. Ambiguous transcript models are not silently assigned a single state. |
| `regulator_abundance` | regulator, population, mean mRNA UMI or missing-prior-only status | 10x mRNA counts do not measure miRNA abundance. |
| `evidence` | accession, source/database, version, assay, cell context, tier, URL | One row per distinct evidence item used by a gene. `not_available` means no externally verified accession. |

The package rejects missing PAS IDs, unknown cell barcodes, cluster mismatches, invalid counts, and duplicate exact cell keys. A PAS with no compatible 3′UTR transcript remains as `no_compatible_3utr_model`. The output `provenance.json` records source paths, sizes, timestamps and SHA-256 where practical. All output CSVs contain only a selected gene to keep the demonstration compact.

Known source discrepancies are detailed in [DATA_INVENTORY.md](DATA_INVENTORY.md). The earlier inventory covers a prior PBMC 1k run; the current MVP takes the completed replica through the `UTRNET_UPSTREAM` configuration variable.
