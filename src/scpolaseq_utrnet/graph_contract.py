"""Framework-independent export contract for stable UTRNet biological graphs.

The contract is deliberately layered on the existing CSV and GraphML outputs.
It does not construct or score graphs and has no PyTorch dependency.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Any
from urllib.parse import quote

import networkx as nx
import pandas as pd

from .metric_registry import METRIC_REGISTRY

GRAPH_CONTRACT_VERSION = "utrnet-graph/v0.1"
BUNDLE_DIRECTORY = "utrnet-graph-v0.1"

REQUIRED_TABLES: dict[str, tuple[str, ...]] = {
    "apa_profiles.csv": (
        "gene_id", "celltype_corrected", "site_id", "umi_sum",
        "total_gene_umi", "pas_proportion", "n_cells",
    ),
    "utr_isoforms.csv": (
        "isoform_id", "gene_id", "site_id", "sequence_length",
        "ambiguity_status", "model_weight",
    ),
    "regulatory_sites.csv": (
        "regulatory_site_id", "isoform_id", "site_id", "regulator",
        "regulator_class", "genomic_start", "genomic_end", "accessibility",
        "evidence_tier",
    ),
    "regulator_abundance.csv": ("population", "regulator", "mean_umi_per_cell"),
    "interaction_scores.csv": (
        "population", "isoform_id", "site_id", "regulatory_site_id",
        "regulator", "regulator_class", "model_weight", "accessibility",
        "evidence_factor", "evidence_tier", "abundance_factor",
        "interaction_score", "score_status",
    ),
}

OPTIONAL_TABLES: dict[str, tuple[str, ...]] = {
    "pas_sites.csv": ("site_id", "start", "end", "strand", "gene_id"),
    "cpa_regulatory_edges.csv": ("source", "factor", "correlation_r", "p_value", "evidence_tier"),
    "cpa_expression.csv": ("population", "factor", "mean_expr", "pct_expressing"),
}

SORT_KEYS: dict[str, tuple[str, ...]] = {
    "apa_profiles.csv": ("gene_id", "celltype_corrected", "site_id"),
    "pas_sites.csv": ("gene_id", "site_id", "start"),
    "utr_isoforms.csv": ("gene_id", "site_id", "isoform_id"),
    "regulatory_sites.csv": ("isoform_id", "site_id", "regulatory_site_id"),
    "regulator_abundance.csv": ("population", "regulator"),
    "interaction_scores.csv": ("population", "isoform_id", "regulatory_site_id", "regulator"),
    "cpa_regulatory_edges.csv": ("source", "target", "factor", "evidence_tier"),
    "cpa_expression.csv": ("population", "factor"),
}

AMBIGUITY_MAP = {
    "unique": "unique",
    "same_sequence_multiple_transcripts": "ambiguous_multiple_models",
    "ambiguous_sequence_models": "ambiguous_multiple_models",
    "ambiguous_multiple_models": "ambiguous_multiple_models",
    "model_dependent": "ambiguous_multiple_models",
    "no_compatible_3utr_model": "no_compatible_3utr_model",
}

RELATIONSHIP_CATEGORIES = {
    "gene_to_pas": "structural",
    "pas_to_isoform": "annotation_derived",
    "isoform_to_regulatory_region": "structural",
    "regulatory_region_to_regulator": "inferred_candidate",
    "regulator_to_gene": "inferred_candidate",
    "population_to_expression": "experimental_observation",
    "population_to_isoform_usage": "experimental_observation",
    "cpa_factor_to_gene": "inferred_candidate",
}
EVIDENCE_CATEGORIES = frozenset({
    "structural", "annotation_derived", "experimental_observation", "inferred_candidate",
})
CORE_BUNDLE_FILES = frozenset((*REQUIRED_TABLES, "apa_metric_registry.json"))


@dataclass(frozen=True)
class GraphValidationReport:
    valid: bool
    contract_version: str
    source_path: str
    bundle_path: str | None
    present_tables: list[str]
    row_counts: dict[str, int]
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    bundle_integrity: str = "not_checked"
    source_consistency: str = "not_checked"
    metric_compatibility: str = "not_checked"


def canonical_entity_id(entity_type: str, source_id: Any) -> str:
    """Build a reversible, delimiter-safe identifier without rewriting source IDs."""
    if entity_type not in {"population", "gene", "pas", "utr", "region", "regulator", "cpa"}:
        raise ValueError(f"unknown graph entity type: {entity_type}")
    value = str(source_id)
    if not value or value.lower() in {"nan", "none"}:
        raise ValueError(f"empty source identifier for {entity_type}")
    return f"{entity_type}:{quote(value, safe='._-')}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_snapshot(directory: Path) -> dict[str, dict[str, Any]]:
    paths = [*directory.glob("*.csv")]
    if (directory / "provenance.json").is_file():
        paths.append(directory / "provenance.json")
    if (directory / "graphs").is_dir():
        paths.extend((directory / "graphs").glob("*.graphml"))
        paths.extend((directory / "graphs").glob("*.json"))
    return {
        path.relative_to(directory).as_posix(): {"sha256": _sha256(path), "bytes": path.stat().st_size}
        for path in sorted(path for path in paths if path.is_file())
    }


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _validate_csv_tables(directory: Path, *, check_manifest: bool) -> tuple[list[str], dict[str, int], list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []
    present: list[str] = []
    counts: dict[str, int] = {}
    frames: dict[str, pd.DataFrame] = {}

    for name, required in REQUIRED_TABLES.items():
        path = directory / name
        if not path.is_file():
            errors.append(f"missing required table: {name}")
            continue
        try:
            frame = _read_csv(path)
        except Exception as exc:  # malformed user-supplied table
            errors.append(f"cannot read {name}: {exc}")
            continue
        present.append(name)
        frames[name] = frame
        counts[name] = len(frame)
        missing = sorted(set(required) - set(frame.columns))
        if missing:
            errors.append(f"{name} missing required columns: {missing}")
        elif name in {"apa_profiles.csv", "utr_isoforms.csv"} and frame.empty:
            errors.append(f"required table has no rows: {name}")

    for name, required in OPTIONAL_TABLES.items():
        path = directory / name
        if not path.is_file():
            continue
        try:
            frame = _read_csv(path)
        except Exception as exc:
            errors.append(f"cannot read optional table {name}: {exc}")
            continue
        present.append(name)
        frames[name] = frame
        counts[name] = len(frame)
        missing = sorted(set(required) - set(frame.columns))
        if missing:
            errors.append(f"optional table {name} missing required columns: {missing}")

    for path in sorted(directory.glob("*.csv")):
        if path.name not in present:
            present.append(path.name)
            try:
                frames[path.name] = _read_csv(path)
                counts[path.name] = len(frames[path.name])
            except Exception as exc:
                errors.append(f"cannot read additional table {path.name}: {exc}")

    required_frames = {name: frames.get(name) for name in REQUIRED_TABLES}
    structurally_complete = all(
        frame is not None and set(REQUIRED_TABLES[name]).issubset(frame.columns)
        for name, frame in required_frames.items()
    )
    if structurally_complete:
        profiles = required_frames["apa_profiles.csv"]
        isoforms = required_frames["utr_isoforms.csv"]
        sites = required_frames["regulatory_sites.csv"]
        abundance = required_frames["regulator_abundance.csv"]
        scores = required_frames["interaction_scores.csv"]

        identity_columns = {
            "apa_profiles.csv": (profiles, ("gene_id", "celltype_corrected", "site_id")),
            "utr_isoforms.csv": (isoforms, ("isoform_id",)),
            "regulatory_sites.csv": (sites, ("regulatory_site_id",)),
            "regulator_abundance.csv": (abundance, ("population", "regulator")),
            "interaction_scores.csv": (scores, ("population", "isoform_id", "regulatory_site_id", "regulator")),
        }
        for name, (frame, key) in identity_columns.items():
            if not set(key).issubset(frame.columns):
                continue
            if frame[list(key)].eq("").any(axis=None):
                errors.append(f"{name} has empty identity fields in {key}")
            if frame.duplicated(list(key)).any():
                errors.append(f"{name} has duplicate identity key {key}")

        allowed_ambiguity = set(AMBIGUITY_MAP) | set(AMBIGUITY_MAP.values())
        if "ambiguity_status" in isoforms:
            unknown = sorted(set(isoforms.ambiguity_status) - allowed_ambiguity)
            if unknown:
                errors.append(f"utr_isoforms.csv has unrecognized ambiguity_status values: {unknown}")

        for source_table, column, target in (
            ("apa_profiles.csv", "site_id", set(isoforms.site_id)),
            ("apa_profiles.csv", "gene_id", set(isoforms.gene_id)),
            ("regulatory_sites.csv", "isoform_id", set(isoforms.isoform_id)),
            ("regulatory_sites.csv", "site_id", set(isoforms.site_id)),
            ("interaction_scores.csv", "regulatory_site_id", set(sites.regulatory_site_id)),
            ("interaction_scores.csv", "isoform_id", set(isoforms.isoform_id)),
        ):
            frame = frames[source_table]
            if column in frame and not set(frame[column]) <= target:
                errors.append(f"{source_table}.{column} contains identifiers absent from its source table")

        isoform_site = set(zip(isoforms.isoform_id, isoforms.site_id))
        for table_name, frame in (("regulatory_sites.csv", sites), ("interaction_scores.csv", scores)):
            if not set(zip(frame.isoform_id, frame.site_id)) <= isoform_site:
                errors.append(f"{table_name} contains an isoform_id/site_id pair not present in utr_isoforms.csv")
        region_tuples = set(zip(sites.regulatory_site_id, sites.isoform_id, sites.site_id, sites.regulator))
        score_tuples = set(zip(scores.regulatory_site_id, scores.isoform_id, scores.site_id, scores.regulator))
        if not score_tuples <= region_tuples:
            errors.append("interaction_scores.csv contains a regulatory_site_id/isoform_id/site_id/regulator tuple absent from regulatory_sites.csv")
        site_gene = set(zip(isoforms.site_id, isoforms.gene_id))
        if not set(zip(profiles.site_id, profiles.gene_id)) <= site_gene:
            errors.append("apa_profiles.csv contains a site_id/gene_id pair absent from utr_isoforms.csv")
        pas_frame = frames.get("pas_sites.csv")
        if pas_frame is not None and {"site_id", "gene_id"}.issubset(pas_frame.columns):
            pas_gene = set(zip(pas_frame.site_id, pas_frame.gene_id))
            if not set(zip(profiles.site_id, profiles.gene_id)) <= pas_gene:
                errors.append("apa_profiles.csv contains a site_id/gene_id pair absent from pas_sites.csv")
        if "gene_id" in scores.columns:
            isoform_gene = set(zip(isoforms.isoform_id, isoforms.gene_id))
            if not set(zip(scores.isoform_id, scores.gene_id)) <= isoform_gene:
                errors.append("interaction_scores.csv contains an isoform_id/gene_id pair absent from utr_isoforms.csv")

        if "evidence_tier" in sites:
            if sites.evidence_tier.eq("").any():
                errors.append("regulatory_sites.csv contains empty evidence_tier values")
        if "score_status" in scores and scores.score_status.eq("").any():
            errors.append("interaction_scores.csv contains empty score_status values")

    if check_manifest:
        manifest_path = directory / "graph_contract.json"
        if not manifest_path.is_file():
            errors.append("missing graph_contract.json")
        else:
            try:
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                if manifest.get("contract_version") != GRAPH_CONTRACT_VERSION:
                    errors.append("graph_contract.json has an unsupported contract_version")
                files = manifest.get("files")
                if not isinstance(files, dict):
                    errors.append("graph_contract.json files must be an object")
                    files = {}
                missing_core = sorted(CORE_BUNDLE_FILES - set(files))
                if missing_core:
                    errors.append(f"graph_contract.json is missing required file hashes: {missing_core}")
                actual_files = {
                    path.relative_to(directory).as_posix()
                    for path in directory.rglob("*") if path.is_file() and path.name != "graph_contract.json"
                }
                if set(files) != actual_files:
                    errors.append("graph_contract.json file inventory does not match bundle contents")
                for name, expected in files.items():
                    if Path(name).is_absolute() or ".." in Path(name).parts or not isinstance(expected, dict) or not isinstance(expected.get("sha256"), str):
                        errors.append(f"manifest has invalid file checksum entry: {name}")
                        continue
                    path = directory / name
                    if not path.is_file():
                        errors.append(f"manifest file is missing: {name}")
                    elif _sha256(path) != expected["sha256"]:
                        errors.append(f"manifest hash mismatch: {name}")
                categories = manifest.get("evidence_categories", [])
                if not isinstance(categories, list) or set(categories) != EVIDENCE_CATEGORIES:
                    errors.append("graph_contract.json evidence_categories do not match the canonical vocabulary")
                relationship_categories = manifest.get("relationship_categories", {})
                if not isinstance(relationship_categories, dict) or not set(relationship_categories.values()) <= EVIDENCE_CATEGORIES:
                    errors.append("graph_contract.json relationship_categories contain an unadvertised category")
            except Exception as exc:
                errors.append(f"cannot validate graph_contract.json: {exc}")
    return sorted(set(present)), counts, errors, warnings


def validate_graph_export(
    directory: str | Path,
    *,
    require_manifest: bool = False,
    exact_directory: bool = False,
    raise_on_error: bool = False,
) -> GraphValidationReport:
    """Validate a legacy export directory or a versioned contract bundle."""
    source = Path(directory)
    bundle = source / BUNDLE_DIRECTORY if not exact_directory and (source / BUNDLE_DIRECTORY).is_dir() else None
    target = bundle or source
    if not target.is_dir():
        report = GraphValidationReport(
            False, GRAPH_CONTRACT_VERSION, str(source), None, [], {},
            [f"graph export directory does not exist: {source}"], [],
        )
        if raise_on_error:
            raise ValueError("\n".join(report.errors))
        return report
    present, counts, errors, warnings = _validate_csv_tables(target, check_manifest=require_manifest or bundle is not None)
    for graphml in sorted((target / "graphs").glob("*.graphml")) if (target / "graphs").is_dir() else []:
        try:
            graph = nx.read_graphml(graphml)
            counts[f"graphs/{graphml.name}.nodes"] = graph.number_of_nodes()
            counts[f"graphs/{graphml.name}.edges"] = graph.number_of_edges()
            for node_id, attrs in graph.nodes(data=True):
                kind = str(attrs.get("kind", ""))
                if not kind:
                    errors.append(f"{graphml.name}: node {node_id!r} has no kind")
            graph_edges = (graph.edges(keys=True, data=True) if graph.is_multigraph()
                           else ((left, right, None, attrs) for left, right, attrs in graph.edges(data=True)))
            for left, right, _, attrs in graph_edges:
                if not attrs.get("kind"):
                    errors.append(f"{graphml.name}: edge ({left!r}, {right!r}) has no kind")
                if not attrs.get("evidence"):
                    errors.append(f"{graphml.name}: edge ({left!r}, {right!r}) has no evidence category")
        except Exception as exc:
            errors.append(f"cannot read GraphML {graphml.name}: {exc}")
    source_consistency = "not_checked"
    metric_compatibility = "not_checked"
    if bundle is not None:
        try:
            manifest = json.loads((bundle / "graph_contract.json").read_text(encoding="utf-8"))
            snapshot = manifest.get("source_snapshot")
            metric_compatibility = str(manifest.get("metric_compatibility", {}).get("status", "not_declared"))
            if isinstance(snapshot, dict):
                current = _source_snapshot(source)
                source_consistency = "current" if current == {
                    name: value.get("sha256") for name, value in snapshot.items()
                } else "stale"
                if source_consistency == "stale":
                    warnings.append("bundle is intact but source_snapshot differs from current source files")
            elif bundle is not None:
                source_consistency = "unavailable"
                warnings.append("bundle integrity is checked; source snapshot is unavailable")
        except (OSError, json.JSONDecodeError, AttributeError):
            source_consistency = "unavailable"
    bundle_integrity = "valid" if (require_manifest or bundle is not None) and not errors else (
        "invalid" if require_manifest or bundle is not None else "not_checked"
    )
    reported_bundle_path = str(bundle or target) if require_manifest and (target / "graph_contract.json").is_file() else (
        str(bundle) if bundle is not None else None
    )
    report = GraphValidationReport(
        not errors, GRAPH_CONTRACT_VERSION, str(source), reported_bundle_path,
        present, counts, sorted(set(errors)), sorted(set(warnings)), bundle_integrity,
        source_consistency, metric_compatibility,
    )
    if report.errors and raise_on_error:
        raise ValueError("Invalid utrnet-graph/v0.1 export:\n- " + "\n- ".join(report.errors))
    return report


def _canonicalize_bundle_table(name: str, frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    if name == "utr_isoforms.csv" and "ambiguity_status" in result:
        result["source_ambiguity_status"] = result["ambiguity_status"]
        result["ambiguity_status"] = result["ambiguity_status"].map(AMBIGUITY_MAP).fillna(result["ambiguity_status"])
    if name == "pas_sites.csv" and "chrom" in result and "chr" not in result:
        result["chr"] = result["chrom"]
    keys = [key for key in SORT_KEYS.get(name, ()) if key in result.columns]
    if keys:
        result = result.sort_values(keys, kind="mergesort", na_position="last").reset_index(drop=True)
    return result


def _entity_inventory(bundle: Path) -> dict[str, list[str]]:
    profiles = _read_csv(bundle / "apa_profiles.csv")
    isoforms = _read_csv(bundle / "utr_isoforms.csv")
    sites = _read_csv(bundle / "regulatory_sites.csv")
    scores = _read_csv(bundle / "interaction_scores.csv")
    entities: dict[str, set[str]] = {
        "population": set(), "gene": set(), "pas": set(), "utr_isoform": set(),
        "regulatory_region": set(), "regulator": set(), "cpa_factor": set(),
    }
    entities["population"].update(canonical_entity_id("population", value) for value in profiles.celltype_corrected)
    entities["gene"].update(canonical_entity_id("gene", value) for value in profiles.gene_id)
    entities["pas"].update(canonical_entity_id("pas", value) for value in profiles.site_id)
    entities["utr_isoform"].update(canonical_entity_id("utr", value) for value in isoforms.isoform_id)
    entities["regulatory_region"].update(canonical_entity_id("region", value) for value in sites.regulatory_site_id)
    entities["regulator"].update(canonical_entity_id("regulator", value) for value in sites.regulator)
    if (bundle / "cpa_expression.csv").is_file():
        cpa = _read_csv(bundle / "cpa_expression.csv")
        if "factor" in cpa.columns:
            entities["cpa_factor"].update(canonical_entity_id("cpa", value) for value in cpa.factor)
    return {kind: sorted(values) for kind, values in sorted(entities.items())}


def _table_state_summary(bundle: Path) -> dict[str, Any]:
    """Summarize raw zeroes separately from evidence-backed zeroes and placeholders."""
    summary: dict[str, Any] = {}
    for path in sorted(bundle.glob("*.csv")):
        frame = _read_csv(path)
        fields: dict[str, Any] = {}
        for column in frame.columns:
            values = frame[column].astype(str)
            missing = values.str.strip().eq("")
            field_state: dict[str, Any] = {"missing": int(missing.sum())}
            numeric = pd.to_numeric(values.where(~missing), errors="coerce")
            if numeric.notna().any():
                raw_zero = numeric == 0
                field_state["raw_numeric_zero"] = int(raw_zero.sum())
                status_column = next((column for column in ("score_status", "coverage_status", "measurement_status")
                                      if column in frame.columns), None)
                if status_column:
                    statuses = frame[status_column].astype(str)
                    placeholder = statuses.isin({"no_observed_gene_umi", "missing_unmeasured", "missing_zero_coverage"})
                    field_state["placeholder_zero"] = int((raw_zero & placeholder).sum())
                    field_state["zero_by_status"] = {
                        str(status): int((raw_zero & (statuses == status)).sum())
                        for status in sorted(statuses[raw_zero].unique())
                    }
                if path.name == "apa_profiles.csv" and column in {"umi_sum", "total_gene_umi", "pas_proportion"}:
                    if column in {"umi_sum", "pas_proportion"}:
                        covered = pd.to_numeric(frame.total_gene_umi, errors="coerce").fillna(0) > 0
                        field_state["placeholder_zero"] = int((raw_zero & ~covered).sum())
                    else:
                        covered = pd.to_numeric(frame.n_cells, errors="coerce").fillna(0) > 0
                    field_state["confirmed_measured_zero"] = int((raw_zero & covered).sum())
            if "ambiguity" in column.lower() or column.endswith("_status"):
                field_state["values"] = {
                    str(key): int(value)
                    for key, value in values.value_counts(dropna=False).sort_index().items()
                }
            fields[column] = field_state
        summary[path.name] = {"rows": len(frame), "fields": fields}
    return summary


def _relationship_counts(bundle: Path) -> dict[str, int]:
    """Count unique source-table relationships represented by the exported tables."""
    profiles = _read_csv(bundle / "apa_profiles.csv")
    isoforms = _read_csv(bundle / "utr_isoforms.csv")
    sites = _read_csv(bundle / "regulatory_sites.csv")
    scores = _read_csv(bundle / "interaction_scores.csv")

    def count(frame: pd.DataFrame, columns: tuple[str, ...]) -> int:
        if not set(columns).issubset(frame.columns):
            return 0
        return int(frame[list(columns)].drop_duplicates().shape[0])

    profile_population_gene = profiles[["celltype_corrected", "gene_id"]].drop_duplicates()
    profile_population_pas = profiles[["celltype_corrected", "site_id"]].drop_duplicates()
    isoform_pas = isoforms[["site_id", "isoform_id"]].drop_duplicates()
    regulator_gene_count = 0
    if len(scores) and {"isoform_id", "regulator"}.issubset(scores.columns):
        # Isoform annotations are the sole source of gene identity for this relation.
        scored_pairs = scores[["isoform_id", "regulator"]].drop_duplicates()
        scored_gene_pairs = scored_pairs.merge(
            isoforms[["isoform_id", "gene_id"]].drop_duplicates(), on="isoform_id", how="inner",
            validate="many_to_many",
        )
        regulator_gene_count = count(scored_gene_pairs, ("regulator", "gene_id"))
    result = {
        "gene_to_pas": count(profiles, ("gene_id", "site_id")),
        "pas_to_isoform": count(isoforms, ("site_id", "isoform_id")),
        "isoform_to_regulatory_region": count(sites, ("isoform_id", "regulatory_site_id")),
        "regulatory_region_to_regulator": count(sites, ("regulatory_site_id", "regulator")),
        "regulator_to_gene": regulator_gene_count,
        "population_to_expression": count(profile_population_gene, ("celltype_corrected", "gene_id")),
        "population_to_isoform_usage": count(profile_population_pas.merge(isoform_pas, on="site_id"),
                                              ("celltype_corrected", "isoform_id")),
    }
    if (bundle / "cpa_regulatory_edges.csv").is_file():
        cpa = _read_csv(bundle / "cpa_regulatory_edges.csv")
        result["cpa_factor_to_gene"] = count(cpa, ("factor", "target"))
    return result


def _detect_repository_commit() -> str | None:
    root = Path(__file__).resolve().parents[2]
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True, timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    value = result.stdout.strip()
    return value or None


def export_graph_contract(
    directory: str | Path,
    *,
    output_directory: str | Path | None = None,
    generation_config: dict[str, Any] | None = None,
    source_repository_commit: str | None = None,
) -> GraphValidationReport:
    """Create a deterministic, self-describing copy of existing graph outputs."""
    source = Path(directory)
    source_report = validate_graph_export(source, exact_directory=True)
    if not source_report.valid:
        raise ValueError("Cannot export invalid graph directory:\n- " + "\n- ".join(source_report.errors))
    initial_source_snapshot = _source_snapshot(source)
    destination = Path(output_directory) if output_directory is not None else source / BUNDLE_DIRECTORY
    source_resolved = source.resolve()
    destination_resolved = destination.resolve()
    if destination_resolved == source_resolved:
        raise ValueError("contract bundle must be a separate directory from legacy graph outputs")
    if destination_resolved in source_resolved.parents:
        raise ValueError("contract bundle must be a separate directory from legacy graph outputs")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or any(part.is_symlink() for part in (destination.parent, *destination.parent.parents)):
        raise ValueError(f"refusing to export through a symlinked destination path: {destination}")
    stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
    backup: Path | None = None
    try:
        for path in sorted(source.glob("*.csv")):
            frame = _canonicalize_bundle_table(path.name, _read_csv(path))
            frame.to_csv(stage / path.name, index=False, lineterminator="\n")
        if (source / "provenance.json").is_file():
            shutil.copyfile(source / "provenance.json", stage / "provenance.json")
        if (source / "graphs").is_dir():
            graph_dir = stage / "graphs"
            graph_dir.mkdir()
            for graphml in sorted((source / "graphs").glob("*.graphml")):
                shutil.copyfile(graphml, graph_dir / graphml.name)
            for json_graph in sorted((source / "graphs").glob("*.json")):
                shutil.copyfile(json_graph, graph_dir / json_graph.name)

        provenance: dict[str, Any] = {}
        provenance_path = stage / "provenance.json"
        if provenance_path.is_file():
            try:
                provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError) as exc:
                raise ValueError(f"invalid provenance.json: {exc}") from exc

        (stage / "apa_metric_registry.json").write_text(
            json.dumps({"registry_version": "apa-metric-registry/v0.1", "metrics": METRIC_REGISTRY},
                       indent=2, sort_keys=True) + "\n", encoding="utf-8",
        )
        file_records = {}
        for path in sorted(p for p in stage.rglob("*") if p.is_file()):
            rel = path.relative_to(stage).as_posix()
            file_records[rel] = {"sha256": _sha256(path), "bytes": path.stat().st_size}

        profiles = _read_csv(stage / "apa_profiles.csv")
        source_snapshot = _source_snapshot(source)
        if source_snapshot != initial_source_snapshot:
            raise ValueError("Source files changed during export; staged publication was cancelled")
        resolved_config = generation_config or provenance.get("generation_config")
        if not resolved_config:
            resolved_config = {
                "status": "unavailable",
                "reason": "legacy provenance did not record the resolved generation configuration",
            }
        manifest = {
        "contract_version": GRAPH_CONTRACT_VERSION,
        "provider": "scpolaseq_utrnet",
        "source_repository": "https://github.com/martings/scpolaseq_utrnet",
        "source_repository_commit": source_repository_commit or _detect_repository_commit(),
        "source_contract": provenance.get("contract", "utrnet-input/v0.1"),
        "source_dataset": provenance.get("source_dataset", "scPolASeq frozen-output graph"),
        "model_version": provenance.get("model", "unknown"),
        "gene_ids": sorted(set(profiles.gene_id)),
        "population_ids": sorted(set(profiles.celltype_corrected)),
        "reference": {
            "genome_fasta": next((str(item.get("path")) for item in provenance.get("sources", [])
                                   if str(item.get("path", "")).endswith((".fa", ".fasta"))), None),
            "annotation_gtf": next((str(item.get("path")) for item in provenance.get("sources", [])
                                     if str(item.get("path", "")).endswith(".gtf")), None),
            "species": provenance.get("species"),
            "genome_build": provenance.get("genome_build"),
        },
        "sample_context": {"sample_id": provenance.get("sample_id")},
        "generation_config": resolved_config,
        "input_files": provenance.get("sources", []),
        "evidence_sources": sorted({
            (str(row.get("database", row.get("source", "unknown"))), str(row.get("version", "unknown")))
            for row in _read_csv(stage / "evidence.csv").to_dict("records")
        }) if (stage / "evidence.csv").is_file() else [],
        "metric_registry": [entry["metric_id"] for entry in METRIC_REGISTRY],
        "entity_id_format": "<entity-type>:<percent-encoded-source-id>; source IDs remain in table columns",
        "entities": _entity_inventory(stage),
        "relationship_categories": RELATIONSHIP_CATEGORIES,
        "relationship_counts": _relationship_counts(stage),
        "evidence_categories": sorted(EVIDENCE_CATEGORIES),
        "metric_compatibility": {"status": "registry_declared", "registry_version": "apa-metric-registry/v0.1"},
        "source_snapshot": source_snapshot,
        "missingness_policy": "raw numeric zeros are retained; zeroes with unavailable coverage are placeholders; confirmed measured zeroes are reported only when coverage supports that interpretation",
        "table_state_summary": _table_state_summary(stage),
        "files": file_records,
        }
        (stage / "graph_contract.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8",
        )
        staged_report = validate_graph_export(stage, require_manifest=True, exact_directory=True)
        if not staged_report.valid:
            raise ValueError("Staged contract bundle failed validation:\n- " + "\n- ".join(staged_report.errors))

        if destination.exists():
            if not destination.is_dir() or destination.is_symlink():
                raise FileExistsError(f"refusing to replace non-directory or symlink destination: {destination}")
            existing_manifest = destination / "graph_contract.json"
            try:
                owned = json.loads(existing_manifest.read_text(encoding="utf-8")).get("contract_version") == GRAPH_CONTRACT_VERSION
            except (OSError, json.JSONDecodeError, AttributeError):
                owned = False
            if not owned:
                raise FileExistsError(f"refusing to replace directory without an owned {GRAPH_CONTRACT_VERSION} manifest")
            backup = destination.parent / f".{destination.name}.backup-{uuid.uuid4().hex}"
            os.replace(destination, backup)
        try:
            os.replace(stage, destination)
        except Exception:
            if backup is not None and backup.exists() and not destination.exists():
                os.replace(backup, destination)
                backup = None
            raise
        # If there was a previous bundle, keep its uniquely named backup. Both directory
        # renames are same-filesystem operations; a crash between them leaves the old
        # complete bundle in the backup path instead of deleting it.
        final_report = validate_graph_export(destination, require_manifest=True, exact_directory=True)
        current_snapshot = _source_snapshot(source)
        source_hashes = {name: value["sha256"] for name, value in source_snapshot.items()}
        source_consistency = "current" if {
            name: value["sha256"] for name, value in current_snapshot.items()
        } == source_hashes else "stale"
        final_warnings = list(final_report.warnings)
        if source_consistency == "stale":
            final_warnings.append("source files changed during export; bundle remains internally validated")
        return replace(final_report, source_path=str(source), bundle_path=str(destination),
                       source_consistency=source_consistency, warnings=sorted(set(final_warnings)))
    finally:
        if stage.exists():
            shutil.rmtree(stage)
