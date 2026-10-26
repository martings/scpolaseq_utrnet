---
name: parameter-sensitivity-planner
description: Plan UTRNet parameter sweeps and identify which outputs or upstream stages a proposed change invalidates.
---

# Parameter sensitivity planner

Input: proposed parameter changes, frozen-output provenance and current analysis metrics.

Separate graph-only, sequence/structural, downstream APA and upstream PAS parameters. Estimate scientific impact, runtime and memory, regenerated files, invalidated outputs and a validation criterion. Start with a dry-run graph-only sweep. Upstream PAS changes may require Sierra or SCAPTURE reruns; planning a rerun does not authorize executing it.

Output: prioritized sweep/backlog and stage invalidation report. Check that each proposed result can be compared against the frozen baseline. Handoff upstream work only after its expected value and authorization are explicit.
