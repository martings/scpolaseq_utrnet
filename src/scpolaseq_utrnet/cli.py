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
from .sensitivity import plan, PARAMETERS, graph_sweep
from .visualization import render_gene
from .graph_contract import export_graph_contract, validate_graph_export


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
        if name in {"ingest", "build-graph", "render", "compare", "sensitivity"}:
            item.add_argument("--gene", help="gene symbol or Ensembl base ID")
        if name == "sensitivity":
            item.add_argument("--dry-run", action="store_true")
            item.add_argument("--change", action="append", choices=sorted(PARAMETERS))

    cerna_item = sub.add_parser("cerna")
    cerna_item.add_argument("--config", required=True)
    cerna_item.add_argument("--top-n", type=int, default=15, help="number of top robust genes")
    cerna_item.add_argument("--genes", nargs="+", help="explicit list of gene symbols")
    cerna_item.add_argument("--output", help="output directory for ceRNA network")

    traj_item = sub.add_parser("trajectory")
    traj_item.add_argument("--config", required=True)
    traj_item.add_argument("--gene", help="gene symbol or Ensembl base ID")
    traj_item.add_argument("--order", nargs="+", required=True, help="ordered list of population states")
    export_item = sub.add_parser("export-contract", help="validate and export a stable utrnet-graph/v0.1 bundle")
    export_item.add_argument("--directory", required=True, help="existing gene output directory")
    export_item.add_argument("--output", help="bundle output directory; defaults inside the gene directory")
    validate_item = sub.add_parser("validate-contract", help="validate an existing graph directory or contract bundle")
    validate_item.add_argument("--directory", required=True)
    args = parser.parse_args(argv)
    if args.command == "export-contract":
        report = export_graph_contract(args.directory, output_directory=args.output)
        print(json.dumps({"valid": report.valid, "contract_version": report.contract_version,
                          "bundle_path": report.bundle_path, "errors": report.errors,
                          "warnings": report.warnings, "bundle_integrity": report.bundle_integrity,
                          "source_consistency": report.source_consistency,
                          "metric_compatibility": report.metric_compatibility}, indent=2))
        return 0 if report.valid else 1
    if args.command == "validate-contract":
        report = validate_graph_export(args.directory, require_manifest=True)
        print(json.dumps({"valid": report.valid, "contract_version": report.contract_version,
                          "bundle_path": report.bundle_path, "present_tables": report.present_tables,
                          "row_counts": report.row_counts, "errors": report.errors,
                          "warnings": report.warnings, "bundle_integrity": report.bundle_integrity,
                          "source_consistency": report.source_consistency,
                          "metric_compatibility": report.metric_compatibility}, indent=2))
        return 0 if report.valid else 1
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
        if args.dry_run:
            print(json.dumps(plan(args.change or list(PARAMETERS)), indent=2))
            return 0
        if args.change:
            parser.error("--change is supported only with --dry-run; the executable sweep uses fixed graph-only scenarios")
        directory = _gene_directory(config)
        scores = pd.read_csv(directory / "interaction_scores.csv")
        target = directory / "sensitivity.csv"
        graph_sweep(scores, config.min_population_cells).to_csv(target, index=False)
        print(target)
        return 0
    if args.command == "cerna":
        from .cerna import export_cerna_analysis
        out_dir = Path(args.output).resolve() if getattr(args, "output", None) else None
        res = export_cerna_analysis(config, genes=getattr(args, "genes", None), top_n=args.top_n, output_dir=out_dir)
        print(json.dumps({k: str(v) for k, v in res.items()}, indent=2))
        return 0
    if args.command == "trajectory":
        from .dynamics import export_trajectory_analysis
        directory = _gene_directory(config)
        res = export_trajectory_analysis(directory, args.order)
        print(json.dumps({k: str(v) for k, v in res.items()}, indent=2))
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
