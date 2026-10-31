# Tahoe-100M source audit (initial metadata pass)

## Fixed source revision and provenance

- Dataset: tahoebio/Tahoe-100M
- Revision inspected: 349e2a55f5538a3e73b32a06c6bdd8e9e87011e8 (HF main at inspection).
- Dataset license declared on the card: CC0-1.0.
- Provenance: Tahoe Bio describes the atlas as produced with Vevo Therapeutics' Mosaic platform; the card cites the 2025 Tahoe-100M bioRxiv preprint.
- No expression or metadata shards were downloaded for this initial pass.

## Observed configuration and storage inventory

The pinned dataset card declares seven Parquet configurations, each with one train split: expression_data, sample_metadata, gene_metadata, drug_metadata, cell_line_metadata, obs_metadata, and pseudobulk_differential_expression. Expression data has 95,624,334 rows; the card reports 1,693,653,078,843 uncompressed bytes and 337,644,770,670 download bytes for that split. The Hugging Face file browser reports 429 GB total repository storage. The pseudobulk DE table is partitioned into 1,026 Parquet files and its directory is reported as 88.9 GB. These figures describe different encodings/aggregations and should not be added as if they were disjoint downloads.

## Observed expression schema

| Field | Type | Interpretation |
| --- | --- | --- |
| genes | sequence<int64> | Sparse per-cell token IDs, aligned to expressions |
| expressions | sequence<float32> | Raw count values; first element is a marker sentinel |
| drug | string | Treatment; DMSO_TF denotes vehicle |
| sample | string | Sample identifier; distinguishes replicate treatments |
| BARCODE_SUB_LIB_ID | string | Unique cell barcode/sub-library key; joins obs_metadata |
| cell_line_id | string | Cellosaurus ID |
| moa-fine | string | Contextual fine-grained MOA annotation |
| canonical_smiles | string | Compound structure |
| pubchem_cid | string | PubChem identifier; empty for DMSO |
| plate | string | Experimental plate |

Map gene tokens through gene_metadata.token_id to ensembl_id and gene_symbol. Discard the first aligned marker pair before mapping or summing. Preserve original token IDs, raw cell keys and raw counts alongside normalized fields. Ensembl IDs are stated as release 109, GRCh38.

## Matching, replication and dose

DMSO_TF controls must match both experimental plate and cell_line_id. Dose is present in sample_metadata.drugname_drugconc as a stringified list of tuples carrying compound, concentration and unit; retain the original field as well as parsed numeric value and unit. The sample identifier is the replicate unit available in the table. Cell counts measure coverage and are not biological replicates. Exact treatment/control sample counts and plate overlap remain to be calculated from pinned metadata and expression indexes.


## Observed metadata table schemas and shard layout

The dataset card documents these metadata columns:

- gene_metadata: gene_symbol, ensembl_id, token_id.
- sample_metadata: sample, plate, mean_gene_count, mean_tscp_count, mean_mread_count, mean_pcnt_mito, drug, drugname_drugconc.
- drug_metadata: drug, targets, moa-broad, moa-fine, human-approved, clinical-trials, gpt-notes-approval, canonical_smiles, pubchem_cid.
- cell_line_metadata: cell_name, Cell_ID_DepMap, Cell_ID_Cellosaur, Organ, Driver_Gene_Symbol, Driver_VarZyg, Driver_VarType, Driver_ProtEffect_or_CdnaEffect, Driver_Mech_InferDM, Driver_GeneType_DM.
- obs_metadata: BARCODE_SUB_LIB_ID is the documented join key; the card does not enumerate its full columns.
- pseudobulk_differential_expression: separate configuration is present; its fields and statistical aggregation are not described in the card and must be captured from the pinned dataset builder.

At the inspected revision the Hugging Face file browser showed 3,388 expression Parquet shards under data/ (338 GB directory total), and the full repository at about 429 GB. The metadata directory was about 91.2 GB because it includes the 88.9 GB pseudobulk DE directory; that DE directory contains 1,026 Parquet shards. The directly listed metadata file sizes were approximately 2.29 GB for obs_metadata, 1.33 MB for gene_metadata, 65.6 kB for sample_metadata, 40.5 kB for drug_metadata and 19 kB for cell_line_metadata. These are browser-reported file sizes, not measured download traffic.

The source README marks the Parquet configs listed above and records the expression train split as 95,624,334 rows. The file browser's shard totals and repo storage are retained as separate storage observations, not interpreted as row counts.

## Other annotations and measurements

- drug_metadata includes targets, moa-broad, moa-fine, approval/trial annotations, SMILES and PubChem IDs. Targets and MOA labels are metadata curated/inferred with GPT-based annotation and compound-source information; they are contextual metadata, not verified direct mechanisms in the assayed cells.
- cell_line_metadata includes Cellosaurus, DepMap, lineage and driver alteration annotations.
- obs_metadata is a separate Parquet table keyed by BARCODE_SUB_LIB_ID. The current card does not expose its complete column schema; runtime audit must record builder feature metadata before use.
- pseudobulk_differential_expression is available as its own configuration, with 1,026 sharded Parquet files. Its current full schema and aggregation definition still require inspection from the pinned builder before benchmark use.
- This is a fixed endpoint perturbation atlas. The source description does not provide time-course measurements or gene-editing Perturb-seq interventions.

## Download budget and next validation

The full expression split is far outside the pilot budget. Do not materialize it or the full pseudobulk DE configuration. The adapter must require an explicit byte cap for cached expression shards and write the selected config, source revision, sample IDs, control matching, bytes fetched, and checksums to the manifest. Metadata-only builder inspection is the first executable audit. Selection remains pending until actual pinned metadata establish two cell lines, 3–5 conditions, matched plate/cell-line controls, QC coverage, and at least two unique sample IDs per treatment/control group.

## Architecture contract mapping

The adapter's output manifest names cell-regulatory-state/v0.1 and preserves source identities, context, perturbation, dose, raw counts, sparse features, and missing modalities. The existing UTRNet ingestion contract is utrnet-input/v0.1 and is not silently substituted. A canonical serialized schema for cell-regulatory-state/v0.1 was not present in the repository at the inspected main revision; field mapping should be reconciled when that shared schema is supplied.
