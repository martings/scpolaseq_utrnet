"""Read-only file inventory for a frozen PBMC run and downstream output."""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path


def inspect_file(path: Path, root: Path, origin: str) -> dict:
    size = path.stat().st_size
    relative = path.relative_to(root).as_posix()
    item = {"origin": origin, "path": f"${{{origin.upper()}}}/{relative}", "relative_path": relative,
            "bytes": size, "sha256": None, "checksum_status": "omitted_large_or_binary",
            "row_count": None, "columns": None, "schema_version": "unversioned", "complete": size > 0}
    if size <= 700_000_000 and path.suffix.lower() in {".tsv", ".csv", ".json", ".bed", ".bedgraph", ".txt", ".h5ad", ".html", ".png", ".log", ".mtx", ".fai", ".gz"}:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
                h.update(block)
        item["sha256"] = h.hexdigest()
        item["checksum_status"] = "calculated"
    if size and path.suffix.lower() in {".tsv", ".csv", ".bed"}:
        with path.open("rb") as stream:
            item["row_count"] = sum(block.count(b"\n") for block in iter(lambda: stream.read(4 * 1024 * 1024), b""))
        if path.suffix.lower() != ".bed":
            item["row_count"] -= 1
            with path.open(newline="") as stream:
                item["columns"] = next(csv.reader(stream, delimiter="\t" if path.suffix.lower() == ".tsv" else ","))
    if size == 0:
        item["validation_findings"] = ["empty_file"]
    return item


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream", type=Path, required=True)
    parser.add_argument("--downstream", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    roots = {"upstream": args.upstream.resolve(), "downstream": args.downstream.resolve()}
    files = [inspect_file(p, root, origin) for origin, root in roots.items()
             for p in sorted(root.rglob("*")) if p.is_file()]
    result = {"inventory_schema": "utrnet-inventory/v0.1",
              "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "source_roots": {k: f"${{{k.upper()}}}" for k in roots}, "files": files,
              "notes": ["TSV/CSV row counts exclude headers; BED counts include every line.",
                        "Checksum omissions are explicit for large and selected binary files.",
                        "Nonzero size is a preliminary completeness flag, not scientific validation."]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"{args.output}: {len(files)} files, {sum(x['sha256'] is not None for x in files)} checksums")


if __name__ == "__main__":
    main()
