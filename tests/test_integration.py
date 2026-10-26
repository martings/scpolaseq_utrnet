from pathlib import Path
import pandas as pd
import pytest

from examples.make_fixture import make_fixture
from scpolaseq_utrnet.config import load_config
from scpolaseq_utrnet.contracts import validate_counts
from scpolaseq_utrnet.ingest import ingest_gene
from scpolaseq_utrnet.ingest import audit_gene_counts
from scpolaseq_utrnet.isoforms import reconstruct_isoforms
from scpolaseq_utrnet.folding import common_region_accessibility, annotate_accessibility
from scpolaseq_utrnet.motifs import scan_mirna_seeds
from scpolaseq_utrnet.graphs import site_retention
from scpolaseq_utrnet.pipeline import build_gene_graphs
from scpolaseq_utrnet.visualization import render_gene


def test_synthetic_end_to_end(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    tables = ingest_gene(config)
    assert audit_gene_counts(config, tables)["exact_key_equal"]
    assert tables["cells"].cell_id.tolist() == ["sample:sample:BC1-1", "sample:sample:BC2-1"]
    assert tables["pas_counts"].umi_count.sum() == 20
    profiles = tables["apa_profiles"].set_index(["celltype_corrected", "site_id"])
    assert profiles.loc[("A", "plus:140"), "pas_proportion"] == pytest.approx(0.2)
    assert profiles.loc[("B", "plus:140"), "pas_proportion"] == pytest.approx(0.8)
    assert profiles.loc[("A", "plus:125"), "usage_ci_low"] <= profiles.loc[("A", "plus:125"), "pas_proportion"]
    assert profiles.loc[("A", "plus:125"), "usage_ci_high"] >= profiles.loc[("A", "plus:125"), "pas_proportion"]
    result = build_gene_graphs(config)
    directory = result["directory"]
    iso = pd.read_csv(directory / "utr_isoforms.csv")
    assert set(iso.ambiguity_status) == {"unique"}
    assert iso.set_index("site_id").loc["plus:125", "sequence_length"] == 15
    assert iso.set_index("site_id").loc["plus:140", "sequence_length"] == 30
    regulatory = pd.read_csv(directory / "regulatory_sites.csv")
    assert set(regulatory.site_id) == {"plus:140"}
    assert set(regulatory.regulator) == {"ELAVL1", "hsa-miR-fixture"}
    retention = pd.read_csv(directory / "site_retention.csv")
    assert set(retention.retention_status) == {"pas_specific"}
    assert set(retention.absent_pas_ids) == {"plus:125"}
    abundance = pd.read_csv(directory / "regulator_abundance.csv")
    assert abundance[abundance.regulator == "hsa-miR-fixture"].mean_umi_per_cell.isna().all()
    scores = pd.read_csv(directory / "interaction_scores.csv")
    assert (scores.score_status == "partial_missing_accessibility_or_abundance").all()
    a = scores[(scores.population == "A") & (scores.regulator == "ELAVL1")].interaction_score.iloc[0]
    b = scores[(scores.population == "B") & (scores.regulator == "ELAVL1")].interaction_score.iloc[0]
    assert b > a
    assert list((directory / "graphs").glob("*.graphml"))
    html, png = render_gene(directory)
    assert html.stat().st_size > 1000 and png.stat().st_size > 1000


def test_minus_strand_and_suspicious_span(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    from dataclasses import replace
    minus = replace(config, gene="MINUS")
    tables = ingest_gene(minus)
    iso, maps = reconstruct_isoforms(tables["pas_sites"], minus.reference_gtf, minus.reference_fasta)
    assert len(iso) == 2
    assert all(len(maps[x.isoform_id]) == x.sequence_length for x in iso.itertuples())
    assert all(maps[x.isoform_id][0] > maps[x.isoform_id][-1] for x in iso.itertuples())
    short, _ = reconstruct_isoforms(tables["pas_sites"], minus.reference_gtf, minus.reference_fasta, max_span=4)
    assert set(short.ambiguity_status) == {"no_compatible_3utr_model"}
    distant = tables["pas_sites"].copy()
    distant.loc[len(distant)] = distant.iloc[0]
    distant.loc[len(distant) - 1, "site_id"] = "minus:distant"
    distant.loc[len(distant) - 1, "pas_reference_id"] = "pas:minus:distant"
    distant.loc[len(distant) - 1, "start"] = 390
    distant.loc[len(distant) - 1, "end"] = 390
    flagged, _ = reconstruct_isoforms(distant, minus.reference_gtf, minus.reference_fasta, max_span=50)
    assert set(flagged.span_status) == {"suspicious_gene_pas_span"}


def test_reject_missing_pas_and_barcode_collision(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    tables = ingest_gene(config)
    invalid = tables["pas_counts"].copy()
    invalid.loc[0, "site_id"] = "missing_site"
    with pytest.raises(ValueError, match="site_id"):
        validate_counts(invalid, tables["cells"], tables["pas_sites"])
    invalid = tables["pas_counts"].copy()
    invalid.loc[0, "cell_barcode"] = "BC1"
    with pytest.raises(ValueError, match="barcode"):
        validate_counts(invalid, tables["cells"], tables["pas_sites"])


def test_grouping_conservation_catches_changed_umi(tmp_path: Path):
    config = load_config(make_fixture(tmp_path))
    tables = ingest_gene(config)
    path = config.upstream / "sierra_quant/cell_type/sample.cell_type.unlabeled.sierra_quant.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.loc[0, "umi_count"] += 1
    frame.to_csv(path, sep="\t", index=False)
    assert audit_gene_counts(config, tables)["status"] == "mismatch"


def test_common_region_accessibility_and_parameter_bounds(tmp_path: Path):
    sites = pd.DataFrame([
        {"regulator": "ELAVL1", "chrom": "chr1", "genomic_start": 10, "genomic_end": 14,
         "strand": "+", "isoform_id": "short", "accessibility": 0.2},
        {"regulator": "ELAVL1", "chrom": "chr1", "genomic_start": 10, "genomic_end": 14,
         "strand": "+", "isoform_id": "long", "accessibility": 0.8},
    ])
    comparison = common_region_accessibility(sites)
    assert comparison.range_accessibility.iloc[0] == pytest.approx(0.6)
    with pytest.raises(ValueError, match="RNAplfold"):
        annotate_accessibility(sites, pd.DataFrame(), tmp_path, window=5, span=4, unpaired=6)


def test_mirna_seed_prior_and_junction_filter(tmp_path: Path):
    fasta = tmp_path / "mature.fa"
    fasta.write_text(">hsa-miR-test\nAGGGGGGGAAAA\n")
    isoforms = pd.DataFrame([{"isoform_id": "iso", "site_id": "pas", "sequence": "AAAAACCCCCCCAAAA"},
                             {"isoform_id": "short", "site_id": "short_pas", "sequence": "AAAAA"}])
    pas = pd.DataFrame([{"site_id": "pas", "chrom": "chr1", "strand": "+"},
                        {"site_id": "short_pas", "chrom": "chr1", "strand": "+"}])
    maps = {"iso": list(range(101, 117)), "short": list(range(101, 106))}
    matches = scan_mirna_seeds(isoforms, maps, pas, fasta)
    assert len(matches) == 1
    assert matches.iloc[0].isoform_id == "iso"
    assert matches.iloc[0].regulator_class == "miRNA"
    assert matches.iloc[0].evidence_tier == "computational_prediction"
    assert (matches.iloc[0].genomic_start, matches.iloc[0].genomic_end) == (106, 112)
    junction = list(range(101, 117))
    junction[9:] = list(range(201, 208))
    assert scan_mirna_seeds(isoforms, {"iso": junction, "short": maps["short"]}, pas, fasta).empty


def test_site_retention_reports_model_ambiguity():
    iso = pd.DataFrame([{"isoform_id": "a", "site_id": "pas1", "sequence_length": 20},
                        {"isoform_id": "b", "site_id": "pas1", "sequence_length": 20},
                        {"isoform_id": "c", "site_id": "pas2", "sequence_length": 10},
                        {"isoform_id": "d", "site_id": "unresolved", "sequence_length": 0}])
    sites = pd.DataFrame([{"isoform_id": "a", "regulator": "RBP", "regulator_class": "ARE-BP",
                           "chrom": "chr1", "genomic_start": 10, "genomic_end": 14, "strand": "+"}])
    row = site_retention(iso, sites).iloc[0]
    assert row.retention_status == "model_dependent"
    assert row.model_dependent_pas_ids == "pas1"
    assert row.absent_pas_ids == "pas2"
    assert row.unresolved_pas_ids == "unresolved"


def test_public_config_requires_explicit_input_paths(monkeypatch):
    monkeypatch.delenv("UTRNET_UPSTREAM", raising=False)
    monkeypatch.delenv("UTRNET_DOWNSTREAM", raising=False)
    config = Path(__file__).resolve().parents[1] / "configs/pbmc1k.yaml"
    with pytest.raises(ValueError, match="unresolved environment variable in upstream"):
        load_config(config)
