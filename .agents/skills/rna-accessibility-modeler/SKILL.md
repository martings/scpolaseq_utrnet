---
name: rna-accessibility-modeler
description: Estimate isoform-aware RNAplfold accessibility and compare common 3′UTR regions across PAS alternatives.
---

# RNA accessibility modeler

Input: validated isoform sequences and site coordinates, RNAplfold executable/version and window, span and unpaired-length parameters.

Cache predictions by sequence and parameter hash. Calculate site accessibility within sequence bounds and compare identical genomic regions across proximal and distal isoforms. Missing executable or boundary probabilities must remain missing; do not fill in a structural claim. Benchmark runtime and memory on a small panel before a broad run.

Output: per-site probabilities, statuses, cache/provenance and common-region differences. Check strand, index, boundary and cache reproducibility. Handoff to the graph builder only after those checks pass.
