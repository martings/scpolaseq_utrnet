import pandas as pd
import pytest

from scpolaseq_utrnet.tahoe_remote import sample_index_from_metadata_rows


def test_metadata_projection_joins_sample_quality_and_counts_cells():
    sample_metadata = pd.DataFrame([{
        "sample": "sample-1",
        "plate": "plate1",
        "drug": "drug-a",
        "drugname_drugconc": "[('drug-a', 0.1, 'uM')]",
        "mean_gene_count": 1200,
        "mean_tscp_count": 5000,
        "mean_mread_count": 10000,
        "mean_pcnt_mito": 4.2,
    }])
    rows = [
        {"sample": "sample-1", "BARCODE_SUB_LIB_ID": "c1", "plate": "plate1", "drug": "drug-a", "cell_line_id": "CVCL_1"},
        {"sample": "sample-1", "BARCODE_SUB_LIB_ID": "c2", "plate": "plate1", "drug": "drug-a", "cell_line_id": "CVCL_1"},
    ]
    result = sample_index_from_metadata_rows(rows, sample_metadata)
    assert len(result) == 1
    assert result.loc[0, "n_cells"] == 2
    assert result.loc[0, "cell_line_id"] == "CVCL_1"
    assert result.loc[0, "drugname_drugconc"] == sample_metadata.loc[0, "drugname_drugconc"]


def test_metadata_projection_rejects_expression_sample_metadata_conflict():
    sample_metadata = pd.DataFrame([{
        "sample": "sample-1", "plate": "plate1", "drug": "drug-a",
        "drugname_drugconc": "[]",
    }])
    row = {"sample": "sample-1", "BARCODE_SUB_LIB_ID": "c1", "plate": "plate2", "drug": "drug-a", "cell_line_id": "CVCL_1"}
    with pytest.raises(ValueError, match="conflict"):
        sample_index_from_metadata_rows([row], sample_metadata)


def test_shared_sample_id_is_grouped_per_cell_line_and_pass_filter():
    sample_metadata = pd.DataFrame([{
        "sample": "shared-sample", "plate": "plate1", "drug": "drug-a",
        "drugname_drugconc": "dose-a", "mean_gene_count": 100,
        "mean_tscp_count": 500, "mean_mread_count": 700, "mean_pcnt_mito": 5,
    }])
    rows = [
        {"sample": "shared-sample", "BARCODE_SUB_LIB_ID": "a1", "plate": "plate1", "drug": "drug-a", "drugname_drugconc": "dose-a", "cell_line": "CVCL_1", "gene_count": 100, "pcnt_mito": 2, "pass_filter": "full"},
        {"sample": "shared-sample", "BARCODE_SUB_LIB_ID": "b1", "plate": "plate1", "drug": "drug-a", "drugname_drugconc": "dose-a", "cell_line": "CVCL_2", "gene_count": 150, "pcnt_mito": 4, "pass_filter": "full"},
        {"sample": "shared-sample", "BARCODE_SUB_LIB_ID": "a2", "plate": "plate1", "drug": "drug-a", "drugname_drugconc": "dose-a", "cell_line": "CVCL_1", "gene_count": 40, "pcnt_mito": 8, "pass_filter": "minimal"},
    ]
    result = sample_index_from_metadata_rows(rows, sample_metadata)
    assert len(result) == 2
    line_a = result[result.cell_line_id == "CVCL_1"].iloc[0]
    assert line_a.n_cells == 1
    assert line_a.n_cells_minimal == 1
    assert line_a.n_cells_obs_total == 2
    assert line_a.mean_gene_count == 100
