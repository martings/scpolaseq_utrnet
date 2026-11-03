# Tahoe adapter PR #1: code ownership and migration inventory

## Source reviewed

Reviewed open draft PR [`martings/scpolaseq_utrnet#1`](https://github.com/martings/scpolaseq_utrnet/pull/1), branch `codex/tahoe-100m-adapter`. The PR is still open and was not merged, rebased, cherry-picked, or changed by this work. No Tahoe files or dataset shards were downloaded during this inventory.

The PR currently adds the following groups of files:

| PR files | What they do | Recommended owner/action |
| --- | --- | --- |
| `src/scpolaseq_utrnet/tahoe.py`, `tahoe_remote.py`, `tahoe_cli.py` | Pinned-source metadata, dose parsing, sample/control selection, projected metadata scan | Move to UTRNNet as a source connector; retain source revision pin, size budget, cell-line/plate controls and explicit metadata-only mode. Keep provider repository focused on biological graph evidence. |
| `tahoe_cache.py`, `tahoe_resources.py` | Sparse cache limits, cache sidecar and memory estimates | Move to UTRNNet's dataset adapter/runtime layer; reuse byte-budget and provenance patterns, then align cache metadata with `cell-regulatory-state/v0.1`. |
| `tahoe_benchmark.py`, `tahoe_gnn.py` | Zero/ridge benchmark and a PyTorch gene-graph predictor | UTRNNet owns model/training/evaluation. Do not add the GNN to the provider: Gemini's feature branch already has a PyG heterograph model, MLP baseline, leakage guards and perturbation-conditioned head. The two GNN abstractions are not interchangeable. |
| `tahoe_smoke.py`, `examples/tahoe_smoke_fixture.json` | Small synthetic smoke path | Migrate only connector/schema tests or fixture details that do not duplicate UTRNNet's model tests. Keep fixture data synthetic. |
| `tests/test_tahoe*.py` | Metadata, cache, baseline and GNN coverage | Split by ownership: source adapter and byte cap tests to UTRNNet's Tahoe connector; model and benchmark tests to UTRNNet experiment tests. No Tahoe-data test belongs in provider CI. |
| `data/tahoe/pilot_selection_manifest.json`, `docs/TAHOE_100M_AUDIT.md`, `docs/CELL_REGULATORY_STATE_V0.1_MAPPING.md` | Metadata-selected pilot and source audit/mapping | Keep a provenance summary only if it informs stable source interpretation; move cohort-selection artifacts and experimental source audit to UTRNNet. Do not treat the provisional mixed full/minimal counts as expression-compatible coverage. |
| `pyproject.toml` Tahoe/GNN extras and `utrnet-tahoe*` scripts | Optional deps and standalone commands | Move optional dependency and entrypoint ownership with the connector/model code. Avoid making PyTorch or Hugging Face libraries provider runtime dependencies. |

## Already implemented in Gemini's read-only feature branch

The branch `feat/m0-m2-utrnnet-hetero-gnn` contains:

- `src/utrnnet/contracts/utrnet_graph_v0_1.py`: validates the five required provider CSVs, optional PAS/CPA tables, GraphML, and conservation into the consumer schema.
- `src/utrnnet/graphs/adapter.py`: converts provider tables into `HeteroData` with stable typed node IDs, missing masks, reverse typed edges and evidence tiers.
- `src/utrnnet/contracts/schema.py`: `cell-regulatory-state/v0.1`, node/edge schemas, evidence tier map and ambiguity codes.
- `tests/gnn/test_upstream_contract_compatibility.py`: builds a synthetic provider output, validates the upstream directory and GraphML, and checks consumer graph counts.
- `docs/architecture/CELL_REGULATORY_STATE_V0_1.md` and `UPSTREAM_INTERFACE_AND_MIGRATION.md`: ownership boundary and staged migration plan.
- Heterogeneous PyG encoder/decoder, trivial non-graph baseline, perturbation-conditioned head and target-leakage checks.

These are read-only observations of Gemini's branch. It is not the default branch, was not checked out locally, and was not modified. The consumer branch should run its integration test against the contract bundle introduced here before considering dependency pinning.

## Reusable provider work

`utrnet-graph/v0.1` required CSV definitions match the five tables in the Gemini validator. This PR adds stricter cross-table referential checks, file checksums, stable entity inventory and a serialized metric registry without importing model libraries. The graph directory contains the biological reference tables needed by UTRNNet but does not claim to supply Tahoe perturbation observations.

## Migration gates and unresolved items

1. Confirm the exact Gemini branch and commit intended for integration; do not target `main` or overwrite its working branch implicitly.
2. Run UTRNNet's existing synthetic `test_upstream_contract_compatibility.py` against a generated contract bundle in a separate UTRNNet-owned clone/environment.
3. Before moving Tahoe connector files, decide whether the source adapter should output `cell-regulatory-state/v0.1` directly or a source-specific schema plus explicit converter. Keep cell counts, sample replicate units, and multi-cell-line sample keys correct.
4. Do not transfer `tahoe_gnn.py` or training code to this provider. The PR's GNN is an experiment stub and should be reconciled with the active UTRNNet model owner.
5. Keep the Tahoe PR draft open until its owner decides whether to close or supersede it; do not merge PR #1 as part of this graph-contract task.
