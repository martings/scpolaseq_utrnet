"""Unit tests for expanded RBP motif scanning and calibrated miRNA seed matching."""
from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest

from scpolaseq_utrnet.contracts import read_tsv
from scpolaseq_utrnet.motifs import scan_motifs, scan_mirna_seeds


def test_rbp_motifs_curated_resource_validity():
    resource = Path(__file__).resolve().parent.parent / "resources/rbp_motifs_curated.tsv"
    assert resource.exists()
    df = read_tsv(resource)
    assert len(df) >= 80
    assert not df.motif_id.duplicated().any()
    required = {"motif_id", "regulator", "regulator_class", "pattern", "evidence_tier", "source", "version", "url"}
    assert required.issubset(df.columns)
    assert set(df.evidence_tier).issubset({
        "motif_heuristic", "computational_prediction", "external_database",
        "experimental_other_context", "experimental_matching_context"
    })


def test_iupac_motif_scanning():
    # Sequence containing:
    # 1. CACRC (SRSF3 motif CNCRC with N=A) at index 5: CACGC
    # 2. UGUAA (CFIm motif UGUAN with N=A) at index 20: UGUAA
    # 3. AUUUA (ARE-BP exact ACGU) at index 35: AUUUA
    seq = "AAAAACACGCAAAAAAAAAAUGUAAAAAAAAAAAAAUUUAAAAAA"
    isoforms = pd.DataFrame([{
        "isoform_id": "iso1",
        "site_id": "pas1",
        "sequence": seq,
    }])
    coord_maps = {"iso1": list(range(1001, 1001 + len(seq)))}
    sites = pd.DataFrame([{
        "site_id": "pas1",
        "chrom": "chr1",
        "strand": "+",
    }])

    motifs_df = pd.DataFrame([
        {"motif_id": "M_SRSF3", "regulator": "SRSF3", "regulator_class": "SR_protein",
         "pattern": "CNCRC", "evidence_tier": "experimental_other_context", "source": "test", "version": "1.0", "url": "test"},
        {"motif_id": "M_CFIM", "regulator": "NUDT21", "regulator_class": "CFIm",
         "pattern": "UGUAN", "evidence_tier": "experimental_other_context", "source": "test", "version": "1.0", "url": "test"},
        {"motif_id": "M_ARE", "regulator": "ELAVL1", "regulator_class": "ARE-BP",
         "pattern": "AUUUA", "evidence_tier": "experimental_other_context", "source": "test", "version": "1.0", "url": "test"},
    ])
    tmp_file = Path("/tmp/test_iupac_motifs.tsv")
    motifs_df.to_csv(tmp_file, sep="\t", index=False)

    matches = scan_motifs(isoforms, coord_maps, sites, tmp_file)
    assert len(matches) == 3
    assert set(matches.regulator) == {"SRSF3", "NUDT21", "ELAVL1"}

    # Verify SRSF3 match at index 5
    srsf3_match = matches[matches.regulator == "SRSF3"].iloc[0]
    assert srsf3_match.isoform_start == 5
    assert srsf3_match.isoform_end == 10
    assert srsf3_match.genomic_start == 1006
    assert srsf3_match.genomic_end == 1010

    # Test invalid alphabet raises ValueError
    bad_motifs = pd.DataFrame([
        {"motif_id": "M_BAD", "regulator": "X", "regulator_class": "Y",
         "pattern": "ACGT123", "evidence_tier": "motif_heuristic", "source": "t", "version": "1", "url": "t"}
    ])
    bad_file = Path("/tmp/test_bad_motifs.tsv")
    bad_motifs.to_csv(bad_file, sep="\t", index=False)
    with pytest.raises(ValueError, match="unsupported motif alphabet"):
        scan_motifs(isoforms, coord_maps, sites, bad_file)


def test_mirna_seed_types(tmp_path: Path):
    fasta = tmp_path / "mature.fa"
    # Mature: 5'-AGGGGGGGAAAA-3'
    # Base 0 = A
    # Bases 1-7 = GGGGGGG -> complement is CCCCCCC
    # 8mer target: CCCCCCCA
    # 7mer-m8 target: CCCCCCC
    fasta.write_text(">hsa-miR-test\nAGGGGGGGAAAA\n")

    # seq_no_8mer: CCCCCCC followed by U -> only 7mer-m8, no 8mer
    seq_no_8mer = "AAAAACCCCCCCUUUUUUUUUU"
    isoforms = pd.DataFrame([{"isoform_id": "iso", "site_id": "pas", "sequence": seq_no_8mer}])
    maps = {"iso": list(range(101, 101 + len(seq_no_8mer)))}
    sites = pd.DataFrame([{"site_id": "pas", "chrom": "chr1", "strand": "+"}])

    matches_7m8 = scan_mirna_seeds(isoforms, maps, sites, fasta, seed_types=("7mer-m8",))
    assert len(matches_7m8) == 1
    assert matches_7m8.iloc[0].matrix_id == "7mer-m8"

    matches_8mer = scan_mirna_seeds(isoforms, maps, sites, fasta, seed_types=("8mer",))
    assert matches_8mer.empty

    # seq_with_8mer: CCCCCCC followed by A -> both 7mer-m8 and 8mer
    seq_with_8mer = "AAAAACCCCCCCAAAAAAAAAA"
    isoforms_8 = pd.DataFrame([{"isoform_id": "iso8", "site_id": "pas", "sequence": seq_with_8mer}])
    maps_8 = {"iso8": list(range(101, 101 + len(seq_with_8mer)))}

    matches_both = scan_mirna_seeds(isoforms_8, maps_8, sites, fasta, seed_types=("7mer-m8", "8mer"))
    assert len(matches_both) == 2
    assert set(matches_both.matrix_id) == {"7mer-m8", "8mer"}
