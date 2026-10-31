"""CLI for pinned Tahoe metadata audits and manifest selection."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from .tahoe import DATASET_REVISION, select_metadata_first_pilot
from .tahoe_remote import inspect_dataset_revision, inspect_file_budget, stream_obs_metadata_sample_index


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="utrnet-tahoe")
    sub = parser.add_subparsers(dest="command", required=True)

    audit = sub.add_parser("audit", help="inspect Hub metadata without downloading shards")
    audit.add_argument("--revision", default=DATASET_REVISION)
    audit.add_argument("--output", type=Path, required=True)
    audit.add_argument("--max-metadata-bytes", type=int, required=True)

    index = sub.add_parser("index", help="stream projected obs_metadata within an explicit byte cap")
    index.add_argument("--output", type=Path, required=True)
    index.add_argument("--max-metadata-bytes", type=int, required=True)

    select = sub.add_parser("select", help="select replicated conditions from a metadata sample index")
    select.add_argument("--sample-index", type=Path, required=True)
    select.add_argument("--output", type=Path, required=True)
    select.add_argument("--cell-lines", type=int, default=2)
    select.add_argument("--conditions", type=int, default=4)
    select.add_argument("--min-treatment-samples", type=int, default=2)
    select.add_argument("--min-control-samples", type=int, default=2)
    select.add_argument("--min-cells", type=int, default=10)

    args = parser.parse_args(argv)
    if args.command == "audit":
        report = inspect_dataset_revision(args.revision)
        report["budget"] = inspect_file_budget(report, max_metadata_bytes=args.max_metadata_bytes)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(args.output)
        return 0 if report["budget"]["within_budget"] else 2

    if args.command == "index":
        frame = stream_obs_metadata_sample_index(max_metadata_bytes=args.max_metadata_bytes)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, index=False)
        summary = {
            "sample_line_rows": len(frame),
            "unique_samples": int(frame["sample"].nunique()),
            "cell_lines": int(frame["cell_line_id"].nunique()),
            "full_expression_compatible_cells": int(frame["n_cells"].sum()),
            "minimal_filter_cells": int(frame["n_cells_minimal"].sum()),
            "obs_metadata_source_bytes": frame.attrs["obs_metadata_source_bytes"],
            "metadata_scan_byte_budget": frame.attrs["metadata_scan_byte_budget"],
            "sample_index": str(args.output),
        }
        print(json.dumps(summary, indent=2))
        return 0

    frame = pd.read_csv(args.sample_index, keep_default_na=False)
    manifest = select_metadata_first_pilot(
        frame,
        cell_lines=args.cell_lines,
        conditions=args.conditions,
        min_treatment_samples=args.min_treatment_samples,
        min_control_samples=args.min_control_samples,
        min_cells=args.min_cells,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
