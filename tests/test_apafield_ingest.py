"""Unit and integration tests for downstream Apafield metrics and priority table ingest."""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.ingest import ingest_gene, ingest_apafield_metrics, ingest_apa_priority
from scpolaseq_utrnet.pipeline import build_gene_graphs


def test_apafield_ingest_with_fixtures(tmp_path: Path):
    config_path = make_fixture(tmp_path)
    config = load_config(config_path)

    # 1. First verify graceful fallback when tables are missing
    apa_df_missing = ingest_apafield_metrics(config, "ENSGPLUS")
    assert apa_df_missing.empty
    priority_df_missing = ingest_apa_priority(config, "ENSGPLUS")
    assert priority_df_missing.empty

    # 2. Create mock downstream tables
    apa_dir = config.downstream / "apa"
    apa_dir.mkdir(parents=True, exist_ok=True)

    si_mag_data = pd.DataFrame([
        {
            "gene_id_base": "ENSGPLUS",
            "celltype": "A",
            "si": 0.2,
            "magnetization_rank": -0.6,
            "magnetization_genomic": -0.5,
            "entropy": 0.72,
        },
        {
            "gene_id_base": "ENSGPLUS",
            "celltype": "B",
            "si": 0.8,
            "magnetization_rank": 0.6,
            "magnetization_genomic": 0.5,
            "entropy": 0.72,
        },
        {
            "gene_id_base": "OTHERGENE",
            "celltype": "A",
            "si": 0.5,
            "magnetization_rank": 0.0,
            "magnetization_genomic": 0.0,
            "entropy": 1.0,
        },
    ])
    si_mag_data.to_csv(config.si_magnetization_table, sep="\t", index=False)

    priority_data = pd.DataFrame([
        {
            "gene_id_base": "ENSGPLUS",
            "group_a": "A",
            "group_b": "B",
            "delta_m_genomic_fixed": -1.0,
            "delta_si": -0.6,
            "endpoint_artifact_suspect": False,
            "is_robust": True,
            "priority_score": 7.5,
            "adj_p_value": 1e-10,
        }
    ])
    priority_data.to_csv(config.apa_priority_table, sep="\t", index=False)

    # 3. Test ingest functions directly
    apa_df = ingest_apafield_metrics(config, "ENSGPLUS")
    assert len(apa_df) == 2
    assert set(apa_df.celltype) == {"A", "B"}
    assert apa_df.loc[apa_df.celltype == "A", "si"].iloc[0] == pytest.approx(0.2)

    priority_df = ingest_apa_priority(config, "ENSGPLUS")
    assert len(priority_df) == 1
    assert priority_df.iloc[0].priority_score == pytest.approx(7.5)

    # 4. Test integration through ingest_gene
    tables = ingest_gene(config)
    profiles = tables["apa_profiles"]
    assert "si" in profiles.columns
    assert "magnetization_genomic" in profiles.columns
    assert "entropy" in profiles.columns

    prof_a = profiles[profiles.celltype_corrected == "A"]
    assert prof_a["si"].iloc[0] == pytest.approx(0.2)
    assert prof_a["magnetization_genomic"].iloc[0] == pytest.approx(-0.5)

    prof_b = profiles[profiles.celltype_corrected == "B"]
    assert prof_b["si"].iloc[0] == pytest.approx(0.8)
    assert prof_b["magnetization_genomic"].iloc[0] == pytest.approx(0.5)

    # 5. Test full pipeline and rewiring table
    result = build_gene_graphs(config)
    directory = result["directory"]
    assert (directory / "apafield_metrics.csv").exists()
    assert (directory / "apa_priority.csv").exists()

    rewiring = pd.read_csv(directory / "rewiring.csv")
    assert not rewiring.empty
    match_edge = rewiring[(rewiring.population_a == "A") & (rewiring.population_b == "B")]
    assert len(match_edge) > 0
    row = match_edge.iloc[0]
    assert row["delta_m_genomic_fixed"] == pytest.approx(-1.0)
    assert row["delta_si"] == pytest.approx(-0.6)
    assert bool(row["endpoint_artifact_suspect"]) is False
    assert bool(row["is_robust"]) is True
    assert row["priority_score"] == pytest.approx(7.5)
