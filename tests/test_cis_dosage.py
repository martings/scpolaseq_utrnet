"""Unit tests for cis-regulatory motif survival field and fractional dosage modeling."""
import pytest
import pandas as pd
import numpy as np

from scpolaseq_utrnet.cis import (
    annotate_motif_survival_probabilities,
    compute_fractional_cis_dosage,
    compute_dynamic_dosage_flux,
)


@pytest.fixture
def synthetic_plus_strand_data():
    # Gene on '+' strand: Stop at 1000, PAS_1 at 1200, PAS_2 at 1500, PAS_3 at 2000
    pas_sites = pd.DataFrame([
        {"site_id": "pas_1", "start": 1200, "end": 1200, "strand": "+"},
        {"site_id": "pas_2", "start": 1500, "end": 1500, "strand": "+"},
        {"site_id": "pas_3", "start": 2000, "end": 2000, "strand": "+"},
    ])

    # Motifs:
    # m1 in cUTR: 1050-1055 (before pas_1)
    # m2 in aUTR_1: 1300-1305 (between pas_1 and pas_2)
    # m3 in aUTR_2: 1700-1705 (between pas_2 and pas_3)
    regulatory_sites = pd.DataFrame([
        {
            "regulatory_site_id": "m1", "gene_id": "GENE1", "regulator": "HuR",
            "regulator_class": "ARE-BP", "motif_id": "ARE1", "chrom": "chr1",
            "genomic_start": 1050, "genomic_end": 1055, "strand": "+"
        },
        {
            "regulatory_site_id": "m2", "gene_id": "GENE1", "regulator": "HuR",
            "regulator_class": "ARE-BP", "motif_id": "ARE2", "chrom": "chr1",
            "genomic_start": 1300, "genomic_end": 1305, "strand": "+"
        },
        {
            "regulatory_site_id": "m3", "gene_id": "GENE1", "regulator": "HuR",
            "regulator_class": "ARE-BP", "motif_id": "ARE3", "chrom": "chr1",
            "genomic_start": 1700, "genomic_end": 1705, "strand": "+"
        },
    ])

    # Two cell states:
    # State A (proximal shift): theta_1 = 0.7, theta_2 = 0.2, theta_3 = 0.1
    # State B (distal shift): theta_1 = 0.1, theta_2 = 0.3, theta_3 = 0.6
    profiles = pd.DataFrame([
        {"site_id": "pas_1", "population": "State_A", "usage": 0.7},
        {"site_id": "pas_2", "population": "State_A", "usage": 0.2},
        {"site_id": "pas_3", "population": "State_A", "usage": 0.1},
        {"site_id": "pas_1", "population": "State_B", "usage": 0.1},
        {"site_id": "pas_2", "population": "State_B", "usage": 0.3},
        {"site_id": "pas_3", "population": "State_B", "usage": 0.6},
    ])

    return pas_sites, regulatory_sites, profiles


def test_plus_strand_survival_field(synthetic_plus_strand_data):
    pas_sites, regulatory_sites, profiles = synthetic_plus_strand_data

    surv_df = annotate_motif_survival_probabilities(
        regulatory_sites, pas_sites, profiles,
        population_col="population", usage_col="usage"
    )

    assert not surv_df.empty
    assert len(surv_df) == 6  # 3 motifs x 2 populations

    # Check State A
    a_m1 = surv_df[(surv_df["population"] == "State_A") & (surv_df["regulatory_site_id"] == "m1")].iloc[0]
    assert a_m1["region_type"] == "cUTR"
    assert a_m1["survival_probability"] == 1.0  # cUTR is constitutive

    a_m2 = surv_df[(surv_df["population"] == "State_A") & (surv_df["regulatory_site_id"] == "m2")].iloc[0]
    assert a_m2["region_type"] == "aUTR"
    # m2 (1305) retained by pas_2 (0.2) + pas_3 (0.1) = 0.3
    assert np.isclose(a_m2["survival_probability"], 0.3)

    a_m3 = surv_df[(surv_df["population"] == "State_A") & (surv_df["regulatory_site_id"] == "m3")].iloc[0]
    assert a_m3["region_type"] == "aUTR"
    # m3 (1705) retained only by pas_3 (0.1) = 0.1
    assert np.isclose(a_m3["survival_probability"], 0.1)

    # Check State B
    b_m2 = surv_df[(surv_df["population"] == "State_B") & (surv_df["regulatory_site_id"] == "m2")].iloc[0]
    # m2 retained by pas_2 (0.3) + pas_3 (0.6) = 0.9
    assert np.isclose(b_m2["survival_probability"], 0.9)


def test_fractional_cis_dosage_and_flux(synthetic_plus_strand_data):
    pas_sites, regulatory_sites, profiles = synthetic_plus_strand_data

    surv_df = annotate_motif_survival_probabilities(
        regulatory_sites, pas_sites, profiles,
        population_col="population", usage_col="usage"
    )

    dosage_df = compute_fractional_cis_dosage(surv_df, population_col="population")
    assert len(dosage_df) == 2  # 1 gene x 1 regulator x 2 populations

    row_a = dosage_df[dosage_df["population"] == "State_A"].iloc[0]
    assert row_a["n_motifs_cutr"] == 1
    # aUTR dosage in A: m2 (0.3) + m3 (0.1) = 0.4
    assert np.isclose(row_a["dosage_autr"], 0.4)
    assert np.isclose(row_a["total_dosage"], 1.4)

    row_b = dosage_df[dosage_df["population"] == "State_B"].iloc[0]
    assert row_b["n_motifs_cutr"] == 1
    # aUTR dosage in B: m2 (0.9) + m3 (0.6) = 1.5
    assert np.isclose(row_b["dosage_autr"], 1.5)
    assert np.isclose(row_b["total_dosage"], 2.5)

    # Flux from A to B: delta = 2.5 - 1.4 = +1.1
    flux_df = compute_dynamic_dosage_flux(
        dosage_df, start_population="State_A", end_population="State_B",
        population_col="population", delta_t=1.0
    )
    assert len(flux_df) == 1
    flux_row = flux_df.iloc[0]
    assert np.isclose(flux_row["delta_dosage"], 1.1)
    assert flux_row["dosage_flux_class"] == "cis_motif_expansion"


def test_minus_strand_survival_field():
    # Gene on '-' strand: transcription runs 3000 -> 1000
    # Stop at 2800, PAS_1 at 2500 (proximal), PAS_2 at 2000, PAS_3 at 1200 (distal)
    pas_sites = pd.DataFrame([
        {"site_id": "pas_1", "start": 2500, "end": 2500, "strand": "-"},
        {"site_id": "pas_2", "start": 2000, "end": 2000, "strand": "-"},
        {"site_id": "pas_3", "start": 1200, "end": 1200, "strand": "-"},
    ])

    regulatory_sites = pd.DataFrame([
        # cUTR: 2600-2605 (before proximal pas_1 at 2500)
        {"regulatory_site_id": "m1", "gene_id": "GENE_MINUS", "regulator": "TTP",
         "regulator_class": "ARE-BP", "motif_id": "TTP1", "chrom": "chr2",
         "genomic_start": 2600, "genomic_end": 2605, "strand": "-"},
        # aUTR: 2200-2205 (between pas_1 at 2500 and pas_2 at 2000)
        {"regulatory_site_id": "m2", "gene_id": "GENE_MINUS", "regulator": "TTP",
         "regulator_class": "ARE-BP", "motif_id": "TTP2", "chrom": "chr2",
         "genomic_start": 2200, "genomic_end": 2205, "strand": "-"},
    ])

    profiles = pd.DataFrame([
        {"site_id": "pas_1", "population": "State_X", "usage": 0.4},
        {"site_id": "pas_2", "population": "State_X", "usage": 0.4},
        {"site_id": "pas_3", "population": "State_X", "usage": 0.2},
    ])

    surv_df = annotate_motif_survival_probabilities(
        regulatory_sites, pas_sites, profiles,
        population_col="population", usage_col="usage"
    )

    m1_row = surv_df[surv_df["regulatory_site_id"] == "m1"].iloc[0]
    assert m1_row["region_type"] == "cUTR"
    assert m1_row["survival_probability"] == 1.0

    m2_row = surv_df[surv_df["regulatory_site_id"] == "m2"].iloc[0]
    assert m2_row["region_type"] == "aUTR"
    # For '-' strand, retained if cleavage <= genomic_start (2200)
    # pas_2 (2000 <= 2200) and pas_3 (1200 <= 2200) satisfy this: 0.4 + 0.2 = 0.6
    assert np.isclose(m2_row["survival_probability"], 0.6)
