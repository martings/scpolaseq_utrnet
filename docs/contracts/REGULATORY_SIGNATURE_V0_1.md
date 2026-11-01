# `regulatory-signature/v0.1`

## Purpose

The regulatory signature is a compact biological evidence bridge. It is intentionally separate from graph topology and from neural-network tensors. It can represent scPolASeq/UTRNet outputs, upstream differential expression or pathway results, Tahoe references, and future bulk, single-cell, spatial, proteomic or chromatin datasets without requiring any one modality.

The formal JSON Schema is [regulatory-signature-v0.1.schema.json](regulatory-signature-v0.1.schema.json); `scpolaseq_utrnet.signature` provides semantic validation and canonical JSON writing.

## Top-level structure

- `dataset`: stable dataset ID and optional repository commit, species and reference genome.
- `context_schema`: open metadata describing fields used by this study.
- `observations`: observation ID, optional biological replicate/sample/library IDs, physical timepoint or other study context. Additional context keys are allowed.
- `modalities`: open modality map with status `measured`, `partially_measured` or `not_measured`, plus optional measurement scale, normalization, and observation unit.
- `measurements`: typed records with stable feature and entity IDs, optional observation, value, unit/scale/normalization, uncertainty, statistical support, source/version, and epistemic/missingness status.

`feature_class` distinguishes gene expression, differential expression, PAS usage, APA metrics, transcript associations, regulatory-site availability, candidate regulatory associations, functional enrichment, pathway activity and perturbation response. These remain different biological/statistical quantities.

## Status and missingness rules

- `measured` with numeric zero is an observed zero.
- `missing_unmeasured` requires `value: null` and means the assay/modality was not measured.
- `missing_zero_coverage` requires `value: null` and `coverage: 0`; it is not a zero proportion.
- `ambiguous` carries unresolved model or metric semantics while preserving the original value when one exists.
- `annotation_derived` and `candidate_inferred` remain separate from `experimental_observation`.
- Functional enrichment records require stable `program_id`, `program_source`, and `program_version`; enrichment also records contributing genes and gene universe size. The signature consumes upstream functional results and does not implement GSEA or Reactome.

Every measurement records its `source`. Statistical estimates can include p-value, adjusted p-value, effect size and uncertainty under `statistical_support` without normalizing their interpretation across methods.

## Fixtures and limits

`examples/regulatory_signature_pbmc.json` is a synthetic PBMC-like context with measured PAS usage, an ambiguous imported SI value, candidate RBP evidence, a pathway result, and explicitly unmeasured small RNA. `examples/regulatory_signature_adipocyte.json` demonstrates a non-immune mouse tissue context with physical timepoints, a measured zero, zero-coverage PAS, and an unresolved transcript model.

This is an interchange schema, not an analysis pipeline. It does not calculate APA metrics, pathway activity, perturbation response, or evidence scores. Species/tissue applicability remains source metadata and should be assessed by each consumer.
