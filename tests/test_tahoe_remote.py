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
