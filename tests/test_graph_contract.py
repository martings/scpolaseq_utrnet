"""Synthetic validation and deterministic serialization for utrnet-graph/v0.1."""
from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import pandas as pd
import pytest

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
    assert set(manifest["relationship_categories"].values()) <= set(manifest["evidence_categories"])
    assert manifest["relationship_categories"]["population_to_expression"] == "experimental_observation"
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
    assert profile_states["umi_sum"]["confirmed_measured_zero"] == int((profiles.umi_sum == 0).sum())
    assert profile_states["pas_proportion"]["confirmed_measured_zero"] == int((profiles.pas_proportion == 0).sum())
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


def test_empty_optional_cpa_table_is_valid_but_headerless_table_is_not(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    columns = ["source", "factor", "correlation_r", "p_value", "evidence_tier"]
    optional = legacy / "cpa_regulatory_edges.csv"
    pd.DataFrame(columns=columns).to_csv(optional, index=False)
    assert validate_graph_export(legacy, exact_directory=True).valid
    optional.write_text("\n")
    report = validate_graph_export(legacy, exact_directory=True)
    assert not report.valid
    assert any("optional table cpa_regulatory_edges.csv" in error for error in report.errors)


def _build_legacy(tmp_path: Path) -> Path:
    return build_gene_graphs(load_config(make_fixture(tmp_path)))["directory"]


def test_missing_required_column_returns_structured_report_without_crashing(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    isoforms_path = legacy / "utr_isoforms.csv"
    isoforms = pd.read_csv(isoforms_path)
    isoforms.drop(columns=["site_id"]).to_csv(isoforms_path, index=False)
    report = validate_graph_export(legacy, exact_directory=True)
    assert not report.valid
    assert any("site_id" in error and "missing required columns" in error for error in report.errors)


def test_compound_score_relationships_are_validated(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    import pandas as pd
    scores_path = legacy / "interaction_scores.csv"
    scores = pd.read_csv(scores_path, keep_default_na=False)
    isoforms = pd.read_csv(legacy / "utr_isoforms.csv", keep_default_na=False)
    sites = pd.read_csv(legacy / "regulatory_sites.csv", keep_default_na=False)
    alternate_site = next(value for value in isoforms.site_id.unique() if value != scores.iloc[0].site_id)
    assert alternate_site in set(isoforms.site_id)
    assert scores.iloc[0].regulatory_site_id in set(sites.regulatory_site_id)
    scores.loc[0, "site_id"] = alternate_site
    scores.to_csv(scores_path, index=False)
    report = validate_graph_export(legacy, exact_directory=True)
    assert not report.valid
    assert any("regulatory_site_id/isoform_id/site_id/regulator tuple" in error for error in report.errors)


@pytest.mark.parametrize(
    ("column", "value", "message"),
    [
        ("population", "unknown_population", "population/site_id pair absent from apa_profiles.csv"),
        ("regulator", "unknown_regulator", "population/regulator pair absent from regulator_abundance.csv"),
    ],
)
def test_score_population_joins_are_validated(tmp_path: Path, column, value, message):
    legacy = _build_legacy(tmp_path)
    scores_path = legacy / "interaction_scores.csv"
    scores = pd.read_csv(scores_path, keep_default_na=False)
    scores.loc[0, column] = value
    scores.to_csv(scores_path, index=False)
    report = validate_graph_export(legacy, exact_directory=True)
    assert not report.valid
    assert any(message in error for error in report.errors)


def test_manifest_cannot_omit_required_file_hashes(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    path = bundle / "graph_contract.json"
    manifest = json.loads(path.read_text())
    manifest["files"].pop("interaction_scores.csv")
    path.write_text(json.dumps(manifest))
    report = validate_graph_export(bundle, require_manifest=True, exact_directory=True)
    assert not report.valid
    assert any("missing required file hashes" in error for error in report.errors)


@pytest.mark.parametrize("field", ["entities", "relationship_counts"])
def test_manifest_entity_metadata_is_required_and_derived_from_tables(tmp_path: Path, field):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    manifest_path = bundle / "graph_contract.json"
    manifest = json.loads(manifest_path.read_text())
    if field == "entities":
        manifest.pop(field)
    else:
        manifest[field]["regulator_to_gene"] = -999
    manifest_path.write_text(json.dumps(manifest))
    report = validate_graph_export(bundle)
    assert not report.valid
    assert any(field in error for error in report.errors)


def test_direct_bundle_path_checks_manifest_hashes_and_metric_registry(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    report = validate_graph_export(bundle)
    assert report.valid, report.errors
    assert report.bundle_integrity == "valid"
    assert report.metric_compatibility == "registry_declared"
    profiles = bundle / "apa_profiles.csv"
    profiles.write_text(profiles.read_text() + "\n")
    damaged = validate_graph_export(bundle)
    assert not damaged.valid
    assert any("hash mismatch: apa_profiles.csv" in error for error in damaged.errors)


def test_cpa_edge_factors_are_in_inventory_without_expression_table(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    (legacy / "cpa_expression.csv").unlink(missing_ok=True)
    pd.DataFrame([{
        "source": "cpa:CPSF1", "target": "gene:ENSGPLUS", "factor": "CPSF1",
        "correlation_r": 0.5, "p_value": 0.01, "evidence_tier": "inferred_candidate",
    }]).to_csv(legacy / "cpa_regulatory_edges.csv", index=False)
    destination = tmp_path / "edge-only-bundle"
    report = export_graph_contract(legacy, output_directory=destination)
    assert report.valid, report.errors
    manifest = json.loads((destination / "graph_contract.json").read_text())
    assert "cpa:CPSF1" in manifest["entities"]["cpa_factor"]


def test_single_population_build_writes_header_complete_rewiring(tmp_path: Path):
    import anndata as ad

    config_path = make_fixture(tmp_path)
    config = load_config(config_path)
    adata = ad.read_h5ad(config.curated_h5ad)
    adata.obs["celltype_corrected"] = "A"
    adata.write_h5ad(config.curated_h5ad)
    built = build_gene_graphs(config)
    scores = pd.read_csv(built["scores"])
    rewiring_header = pd.read_csv(built["directory"] / "rewiring.csv", nrows=0).columns
    assert scores.population.nunique() == 1
    assert {
        "population_a", "population_b", "regulator", "genomic_region", "score_a", "score_b",
        "delta_score", "delta_ci_low", "delta_ci_high", "availability_change",
    } <= set(rewiring_header)


def test_stale_valid_bundle_does_not_validate_modified_source(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    source_table = legacy / "apa_profiles.csv"
    source_table.write_text("gene_id\nmalformed\n")
    report = validate_graph_export(legacy, exact_directory=True)
    assert not report.valid
    assert any("missing required columns" in error for error in report.errors)
    assert (bundle / "graph_contract.json").is_file()


def test_reexport_rebuilds_corrupt_owned_bundle_and_retains_backup(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    (bundle / "apa_profiles.csv").write_text("corrupt")
    report = export_graph_contract(legacy)
    assert report.valid, report.errors
    assert validate_graph_export(bundle, require_manifest=True, exact_directory=True).valid
    backups = list(legacy.glob(f".{BUNDLE_DIRECTORY}.backup-*"))
    assert len(backups) == 1
    assert (backups[0] / "graph_contract.json").is_file()


def test_failed_export_preserves_existing_bundle_byte_for_byte(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    bundle = legacy / BUNDLE_DIRECTORY
    old_manifest = (bundle / "graph_contract.json").read_bytes()
    (legacy / "apa_profiles.csv").write_text("bad\n")
    import pytest
    with pytest.raises(ValueError, match="Cannot export invalid"):
        export_graph_contract(legacy)
    assert (bundle / "graph_contract.json").read_bytes() == old_manifest


def test_source_freshness_is_reported_separately_from_bundle_integrity(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    report = validate_graph_export(legacy)
    assert report.valid and report.bundle_integrity == "valid"
    assert report.source_consistency == "current"
    (legacy / "apa_profiles.csv").write_text((legacy / "apa_profiles.csv").read_text() + "\n")
    stale = validate_graph_export(legacy)
    assert stale.valid and stale.bundle_integrity == "valid"
    assert stale.source_consistency == "stale"


def test_reexport_recovers_complete_resolved_generation_config(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    provenance = json.loads((legacy / "provenance.json").read_text())
    manifest = json.loads((legacy / BUNDLE_DIRECTORY / "graph_contract.json").read_text())
    assert manifest["generation_config"] == provenance["generation_config"]
    assert {"min_pas_umi", "fold_window", "fold_span", "fold_unpaired"} <= set(manifest["generation_config"])
    rebuilt = tmp_path / "re-exported"
    assert export_graph_contract(legacy, output_directory=rebuilt).valid
    rebuilt_manifest = json.loads((rebuilt / "graph_contract.json").read_text())
    assert rebuilt_manifest["generation_config"] == manifest["generation_config"]


def test_export_does_not_clean_unrelated_interrupted_staging_directory(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    interrupted = legacy / f".{BUNDLE_DIRECTORY}.staging-interrupted"
    interrupted.mkdir()
    marker = interrupted / "preserve-me.txt"
    marker.write_text("staging may belong to another interrupted process")
    assert export_graph_contract(legacy).valid
    assert marker.read_text() == "staging may belong to another interrupted process"


def test_score_zero_summary_keeps_placeholder_status_distinct(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    import pandas as pd
    scores_path = legacy / "interaction_scores.csv"
    scores = pd.read_csv(scores_path, keep_default_na=False)
    assert not scores.empty
    scores.loc[0, "interaction_score"] = 0
    scores.loc[0, "score_status"] = "no_observed_gene_umi"
    scores.to_csv(scores_path, index=False)
    report = export_graph_contract(legacy, output_directory=tmp_path / "bundle")
    assert report.valid
    manifest = json.loads((tmp_path / "bundle" / "graph_contract.json").read_text())
    state = manifest["table_state_summary"]["interaction_scores.csv"]["fields"]["interaction_score"]
    assert state["raw_numeric_zero"] >= 1
    assert state["placeholder_zero"] >= 1
    assert "confirmed_measured_zero" not in state


def test_regulator_gene_relationship_count_uses_explicit_isoform_gene_join(tmp_path: Path):
    legacy = _build_legacy(tmp_path)
    manifest = json.loads((legacy / BUNDLE_DIRECTORY / "graph_contract.json").read_text())
    assert manifest["relationship_counts"]["regulator_to_gene"] > 0
