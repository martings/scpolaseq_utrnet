"""Minimal UTRNet command-line interface."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import pandas as pd

from .config import load_config
from .ingest import ingest_gene, write_ingestion
from .pipeline import build_gene_graphs
from .graphs import compare_populations
from .sensitivity import plan, PARAMETERS
from .visualization import render_gene


def _gene_directory(config) -> Path:
    from .ingest import _gene_id, load_cells
    cells, adata = load_cells(config)
    try:
        gene_id, _ = _gene_id(config, adata)
    finally:
        adata.file.close()
    return config.output / gene_id


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="utrnet")
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ["audit", "ingest", "candidates", "build-graph", "compare", "render", "sensitivity"]:
        item = sub.add_parser(name)
        item.add_argument("--config", required=True)
        if name in {"ingest", "build-graph", "render", "compare"}:
            item.add_argument("--gene", help="gene symbol or Ensembl base ID")
        if name == "sensitivity":
            item.add_argument("--dry-run", action="store_true")
            item.add_argument("--change", action="append", choices=sorted(PARAMETERS))
    args = parser.parse_args(argv)
    config = load_config(args.config)
    if getattr(args, "gene", None):
        from dataclasses import replace
        config = replace(config, gene=args.gene)

    if args.command == "audit":
        required = [config.curated_h5ad, config.reference_gtf, config.reference_fasta,
                    config.upstream / "pas_reference/pas_reference.tsv",
                    config.upstream / "labels/cell_annotations.tsv",
                    config.upstream / "sierra_quant/cluster"]
        status = {str(p): {"exists": p.exists(), "nonempty": p.exists() and (p.is_dir() or p.stat().st_size > 0)} for p in required}
        report = {"contract": "utrnet-input/v0.1", "sources": status}
        if all(x["nonempty"] for x in status.values()):
            from .ingest import audit_gene_counts
            tables = ingest_gene(config)
            report["gene"] = str(tables["pas_sites"].gene_name.iloc[0])
            report["cells_annotated"] = len(tables["cells"])
            report["cells_curated"] = int(tables["cells"].curated.sum())
            report["pas_sites"] = len(tables["pas_sites"])
            report["sierra_grouping_conservation"] = audit_gene_counts(config, tables)
        print(json.dumps(report, indent=2))
        conservation = report.get("sierra_grouping_conservation", {}).get("status")
        return 0 if all(x["nonempty"] for x in status.values()) and conservation in {"pass", "missing_cell_type_grouping"} else 1
    if args.command == "candidates":
        table = pd.read_csv(config.downstream / "apa/celltype_specific_APA_priority_table.tsv", sep="\t", keep_default_na=False)
        selected = table[(table.priority_class == "robust") & (table.gene_id_base != "NA") &
                         (table.basis_span_bp <= config.max_utr_span)]
        result = selected.sort_values("priority_score", ascending=False)[
            ["gene_id_base", "gene_name", "comparison", "basis_span_bp", "priority_score"]].head(25)
        print(result.to_string(index=False))
        return 0
    if args.command == "ingest":
        tables = ingest_gene(config)
        directory = _gene_directory(config)
        write_ingestion(tables, directory)
        print(directory)
        return 0
    if args.command == "build-graph":
        print(json.dumps({k: str(v) for k, v in build_gene_graphs(config).items()}, indent=2))
        return 0
    if args.command == "compare":
        directory = _gene_directory(config)
        scores = pd.read_csv(directory / "interaction_scores.csv")
        result = compare_populations(scores)
        result.to_csv(directory / "rewiring.csv", index=False)
        print(directory / "rewiring.csv")
        return 0
    if args.command == "render":
        print("\n".join(str(x) for x in render_gene(_gene_directory(config))))
        return 0
    if args.command == "sensitivity":
        if not args.dry_run:
            parser.error("sensitivity currently requires --dry-run")
        print(json.dumps(plan(args.change or list(PARAMETERS)), indent=2))
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
