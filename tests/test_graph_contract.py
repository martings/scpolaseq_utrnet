"""Synthetic validation and deterministic serialization for utrnet-graph/v0.1."""
from __future__ import annotations

import json
from pathlib import Path

import networkx as nx

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.graph_contract import (
    BUNDLE_DIRECTORY,
    canonical_entity_id,
    export_graph_contract,
    validate_graph_export,
)
from scpolaseq_utrnet.pipeline import build_gene_graphs


def _graph_identity(path: Path):
    graph = nx.read_graphml(path)
    nodes = sorted((str(node), str(attrs.get("kind", ""))) for node, attrs in graph.nodes(data=True))
    edges = sorted((str(left), str(right), str(attrs.get("kind", "")), str(attrs.get("evidence", "")))
                   for left, right, attrs in graph.edges(data=True))
    return nodes, edges


def test_synthetic_graph_bundle_validates_and_preserves_relations(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    built = build_gene_graphs(config)
    legacy = built["directory"]
    bundle = legacy / BUNDLE_DIRECTORY

    report = validate_graph_export(bundle, require_manifest=True)
    assert report.valid, report.errors
    assert report.contract_version == "utrnet-graph/v0.1"
    assert set(report.present_tables) >= {
        "apa_profiles.csv", "pas_sites.csv", "utr_isoforms.csv", "regulatory_sites.csv",
        "regulator_abundance.csv", "interaction_scores.csv", "evidence.csv",
    }

    manifest = json.loads((bundle / "graph_contract.json").read_text())
    assert manifest["relationship_categories"]["pas_to_isoform"] == "annotation_derived"
    assert manifest["relationship_categories"]["regulatory_region_to_regulator"] == "inferred_candidate"
    assert manifest["relationship_counts"]["gene_to_pas"] > 0
    assert manifest["relationship_counts"]["pas_to_isoform"] > 0
    assert manifest["table_state_summary"]["apa_profiles.csv"]["fields"]["pas_proportion"]["missing"] == 0
    assert manifest["entities"]["gene"] == ["gene:ENSGPLUS"]
    assert manifest["entities"]["pas"]
    assert manifest["entities"]["utr_isoform"]

    legacy_graphs = sorted((legacy / "graphs").glob("*.graphml"))
    bundle_graphs = sorted((bundle / "graphs").glob("*.graphml"))
    assert [p.name for p in legacy_graphs] == [p.name for p in bundle_graphs]
    for source, exported in zip(legacy_graphs, bundle_graphs):
        assert _graph_identity(source) == _graph_identity(exported)


def test_contract_export_is_repeatable_and_keeps_raw_ambiguity(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    built = build_gene_graphs(config)
    legacy = built["directory"]
    second = tmp_path / "second-bundle"
    first = legacy / BUNDLE_DIRECTORY
    first_manifest = json.loads((first / "graph_contract.json").read_text())
    report = export_graph_contract(
        legacy,
        output_directory=second,
        generation_config=first_manifest["generation_config"],
    )
    assert report.valid, report.errors

    for rel in ("apa_profiles.csv", "pas_sites.csv", "utr_isoforms.csv", "regulatory_sites.csv",
                "regulator_abundance.csv", "interaction_scores.csv", "graph_contract.json",
                "apa_metric_registry.json"):
        assert (first / rel).read_bytes() == (second / rel).read_bytes()

    import pandas as pd
    source_isoforms = pd.read_csv(legacy / "utr_isoforms.csv", keep_default_na=False)
    exported_isoforms = pd.read_csv(first / "utr_isoforms.csv", keep_default_na=False)
    assert "source_ambiguity_status" in exported_isoforms
    assert set(exported_isoforms.source_ambiguity_status) == set(source_isoforms.ambiguity_status)
    assert set(exported_isoforms.ambiguity_status) <= {
        "unique", "ambiguous_multiple_models", "no_compatible_3utr_model",
    }


def test_contract_summary_distinguishes_measured_zero_from_missing(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    built = build_gene_graphs(config)
    bundle = built["directory"] / BUNDLE_DIRECTORY
    manifest = json.loads((bundle / "graph_contract.json").read_text())
    import pandas as pd
    profiles = pd.read_csv(bundle / "apa_profiles.csv")
    profile_states = manifest["table_state_summary"]["apa_profiles.csv"]["fields"]
    assert profile_states["umi_sum"]["measured_zero"] == int((profiles.umi_sum == 0).sum())
    assert profile_states["pas_proportion"]["measured_zero"] == int((profiles.pas_proportion == 0).sum())
    isoform_states = manifest["table_state_summary"]["utr_isoforms.csv"]["fields"]
    assert "values" in isoform_states["ambiguity_status"]
    assert sum(isoform_states["ambiguity_status"]["values"].values()) > 0


def test_canonical_ids_are_stable_and_delimiter_safe():
    assert canonical_entity_id("gene", "ENSG:1") == "gene:ENSG%3A1"
    assert canonical_entity_id("gene", "ENSG:1") == canonical_entity_id("gene", "ENSG:1")


def test_validator_reports_missing_required_table(tmp_path: Path):
    report = validate_graph_export(tmp_path)
    assert not report.valid
    assert any("apa_profiles.csv" in error for error in report.errors)
