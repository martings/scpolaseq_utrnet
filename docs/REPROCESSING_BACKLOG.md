# Reprocessing and sensitivity backlog

The first executable sweep is graph-only and never reruns STARsolo, Sierra or SCAPTURE. After `utrnet build-graph`, run `utrnet sensitivity --config configs/pbmc1k.yaml --gene SAT1` to write `sensitivity.csv` from frozen interaction factors. It compares the baseline with half/double miRNA evidence weight, a conservative 0.5 factor for missing accessibility, and a stricter minimum cell count. It reports candidate-edge score sums, rank correlation, top-ten overlap, and population inclusion. These scenarios test model assumptions; they are not biological replicates or calibrated significance tests. `utrnet sensitivity --dry-run` separately reports invalidated outputs before any parameter change.

| Parameter group | Examples | Scientific impact | Runtime/memory | Upstream regeneration | Validation before accepting |
| --- | --- | --- | --- | --- | --- |
| Graph only | evidence/score weights, display threshold | Edge ranking and visible rewiring | seconds, small | no | ranking stability and per-factor review |
| Sequence/structure | motif library/version, RNAplfold window/span/unpaired length | site membership and accessibility | minutes to hours; sequence cache | no | sequence-coordinate checks, null motif enrichment, cache reproducibility |
| Downstream APA | min population cells, effective UMI, PAS proportion estimator | usage uncertainty and sparse-population inclusion | minutes; table-scale | no | count conservation, coverage, bootstrap sensitivity |
| Upstream PAS | merge distance, catalog source, SCAPTURE activation | site identities and all dependent APA results | hours to days; large | yes | predeclared gain in supported PAS and external/annotation concordance |

Priorities: (1) exact Sierra grouping equivalence and PAS catalog duplicate policy; (2) transcript/PAS validation for IFITM3 and ZSWIM6; (3) experimental evidence resources and licenses; (4) RNAplfold benchmark; (5) null models and replicate-aware uncertainty; (6) only then decide whether upstream PAS settings merit a rerun. No upstream rerun is authorized by this backlog.
