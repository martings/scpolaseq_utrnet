# scPolASeq-UTRNet

`scpolaseq_utrnet` builds exploratory, population-specific 3′UTR regulatory graphs from frozen scPolASeq PAS counts. The PBMC 1k MVP uses the completed `pbmc_1k_v3` replica and its technically concordant downstream notebook output. It does not rerun alignment, Sierra, SCAPTURE, or PAS discovery.

## PBMC 1k quick start on Deepthought

```bash
python -m pip install -e '.[test]'
export UTRNET_UPSTREAM=/path/to/completed/scPolASeq/pbmc1k
export UTRNET_DOWNSTREAM=/path/to/executed/downstream/outputs/pbmc1k
utrnet audit --config configs/pbmc1k.yaml
utrnet candidates --config configs/pbmc1k.yaml
utrnet ingest --config configs/pbmc1k.yaml --gene SAT1
utrnet build-graph --config configs/pbmc1k.yaml --gene SAT1
utrnet compare --config configs/pbmc1k.yaml --gene SAT1
utrnet render --config configs/pbmc1k.yaml --gene SAT1
utrnet build-graph --config configs/pbmc1k.yaml --gene CORO1B
utrnet render --config configs/pbmc1k.yaml --gene CORO1B
utrnet sensitivity --config configs/pbmc1k.yaml --dry-run
utrnet sensitivity --config configs/pbmc1k.yaml --gene SAT1
```

Set the two input paths to the completed scPolASeq run and its executed downstream notebook output before running the commands. The configuration rejects unset path variables. Outputs are written under `outputs/pbmc1k/<gene ID>/` and excluded from Git. The HTML explorer, static PNG, normalized CSV tables, Cytoscape JSON, GraphML, and provenance are generated locally. See [PBMC validation](docs/PBMC1K_VALIDATION.md).

SAT1 and CORO1B form the first two-gene panel because both have observed alternative PAS and compatible 3′UTR sequence models. CORO1B also has ambiguous alternative sequence models, which stay flagged. `site_retention.csv` reports sites shared across compatible PAS, restricted to some PAS, or dependent on the selected transcript model. The bundled RBP motif subset is reproducibly extracted from ATtRACT with a pinned archive checksum. The bundled MirGeneDB human mature sequences supply exact miRNA seed-match predictions. IFITM3 remains an APA observation but its highlighted proximal PAS lacks a compatible annotated 3′UTR model in this audit; it is not used to claim a regulatory switch. ZSWIM6 remains excluded pending its large-span mapping review.

## Run without Deepthought

```bash
python examples/make_fixture.py /tmp/utrnet-fixture
utrnet build-graph --config /tmp/utrnet-fixture/config.yaml
utrnet render --config /tmp/utrnet-fixture/config.yaml
python -m pytest -q
```

The fixture contains two cells, two strands, exact barcodes with `-1` suffixes, a synthetic 3′UTR motif retained only in a distal isoform, and a prior-only miRNA site. It is synthetic and makes no biological claim.

## Evidence and interpretation

Observed: Sierra UMI counts, curated cell identities, and mRNA expression. Annotation-derived: candidate transcript/PAS associations and reconstructed 3′UTR sequences. External experimental evidence: ATtRACT's ELAVL1 and ZFP36 AUUUA motif records from other contexts. Predicted: matches of those motifs and MirGeneDB-derived exact miRNA seeds to SAT1 isoforms, plus interaction scores. Missing by default: RNAplfold accessibility if the executable is absent, measured RBP/AGO occupancy at SAT1 in PBMC, and endogenous miRNA abundance. The optional miRNA genomic-site import uses a versioned TSV; both miRNA evidence modes remain prior-only without small-RNA measurements.

Scores expose usage, model weight, motif presence, accessibility, evidence factor, and abundance factor separately. Missing accessibility or abundance yields an explicitly labeled partial score. `population_metrics.csv` summarizes candidate-network size and weight concentration. The executable sensitivity command writes `sensitivity.csv` with five graph-only scenarios; `--dry-run` still reports which stages a parameter change would invalidate. The explorer displays up to eight representative regulator patterns for legibility; its site-retention and interaction tables and exported graphs retain the full candidate set. Results compare lineages within one biological sampling context; cells are not independent biological replicates.

## Repository layout

- `src/scpolaseq_utrnet/`: typed ingestion, isoform, evidence, folding, graph and visualization modules.
- `configs/pbmc1k.yaml`: reproducible local input configuration.
- `resources/are_motifs.tsv`: two checksum-pinned ATtRACT records with PubMed, assay and context columns; [resource notes](docs/REGULATORY_RESOURCES.md).
- `resources/mirgenedb3_hsa_mature.fa`: checksum-pinned human mature miRNA sequences for seed-only site predictions.
- `examples/make_fixture.py`, `tests/`: synthetic demonstration and integration checks.
- `docs/`: prior-run inventory, data contract, scientific model, current PBMC validation, and reprocessing backlog.
- `.agents/skills/`: seven task-specific scientific and engineering skills.

No large datasets, reference FASTA/GTF, BAMs, private configuration or credentials belong in Git.

## Tahoe-100M perturbation adapter

The separate Tahoe adapter pins the source revision and keeps all network access explicit. Install its optional tools and inspect metadata without fetching expression shards. The checked-in smoke fixture is synthetic and can be run without network access:

```bash
python -m pip install -e '.[test,tahoe,gnn]'
utrnet-tahoe-smoke
utrnet-tahoe audit --revision 349e2a55f5538a3e73b32a06c6bdd8e9e87011e8 --max-metadata-bytes 2500000000 --output outputs/tahoe/source-audit.json
```

Build a sample index by streaming only projected obs_metadata fields, with a byte cap checked before the scan, then choose the pilot from observed sample, cell-line, dose, QC and plate-matched DMSO metadata:

```bash
utrnet-tahoe index --max-metadata-bytes 2500000000 --output outputs/tahoe/sample-index.csv
utrnet-tahoe select --sample-index outputs/tahoe/sample-index.csv --conditions 4 --output outputs/tahoe/selection-manifest.json
```

The selector uses unique sample IDs as replicate units and cell counts only as coverage. The committed candidate manifest was selected from obs_metadata; its coverage totals include both `full` and `minimal` filter categories. Validate full-only sample-line counts against expression before caching or training. Use tahoe_cache.write_sparse_cache with a positive max_bytes value to save a selected CSR pilot matrix and its identity sidecar. It does not treat Tahoe as time-course or gene-editing Perturb-seq data. See [the source audit](docs/TAHOE_100M_AUDIT.md) and [the contract mapping](docs/CELL_REGULATORY_STATE_V0.1_MAPPING.md). No expression or pseudobulk shards are downloaded by the audit command.
