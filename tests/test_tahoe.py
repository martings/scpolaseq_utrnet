import numpy as np
import pandas as pd
import pytest

from scpolaseq_utrnet.tahoe import (
    expression_row_to_pairs,
    make_gene_map,
    parse_concentrations,
    select_metadata_first_pilot,
    sparse_counts,
    validate_group_partition,
)


def test_token_mapping_and_marker_exclusion_preserve_counts():
    genes = pd.DataFrame({
        "token_id": [10, 11], "gene_symbol": ["A", "B"],
        "ensembl_id": ["ENSG000001", "ENSG000002"],
    })
    row = {"genes": [0, 10, 11], "expressions": [-2, 3, 7]}
    assert expression_row_to_pairs(row, make_gene_map(genes), marker_token_id=0) == [
        ("ENSG000001", 3), ("ENSG000002", 7),
    ]
    matrix, cells, gene_ids = sparse_counts(
        [{**row, "BARCODE_SUB_LIB_ID": "cell-1"}], genes, marker_token_id=0
    )
    assert cells == ["cell-1"]
    assert gene_ids == ["ENSG000001", "ENSG000002"]
    assert matrix.sum() == 10
    assert matrix.toarray().tolist() == [[3, 7]]


def test_expression_validation_rejects_unknown_duplicate_and_noninteger_tokens():
    genes = pd.DataFrame({
        "token_id": [10], "gene_symbol": ["A"], "ensembl_id": ["ENSG000001"],
    })
    mapping = make_gene_map(genes)
    with pytest.raises(ValueError, match="unmapped"):
        expression_row_to_pairs({"genes": [0, 11], "expressions": [-2, 1]}, mapping)
    with pytest.raises(ValueError, match="duplicate"):
        expression_row_to_pairs({"genes": [0, 10, 10], "expressions": [-2, 1, 2]}, mapping)
    with pytest.raises(ValueError, match="invalid raw count"):
        expression_row_to_pairs({"genes": [0, 10], "expressions": [-2, 1.5]}, mapping)


def test_duplicate_cell_id_and_bad_marker_are_rejected():
    genes = pd.DataFrame({
        "token_id": [10], "gene_symbol": ["A"], "ensembl_id": ["ENSG000001"],
    })
    row = {"genes": [0, 10], "expressions": [-2, 1], "BARCODE_SUB_LIB_ID": "x"}
    with pytest.raises(ValueError, match="duplicate BARCODE"):
        sparse_counts([row, row], genes)
    with pytest.raises(ValueError, match="marker"):
        expression_row_to_pairs({"genes": [5, 10], "expressions": [-2, 1]}, make_gene_map(genes), marker_token_id=0)


def test_concentration_parser_retains_value_and_unit():
    parsed = parse_concentrations("[('compound', 0.05, 'uM')]")
    assert parsed == [{"drug": "compound", "concentration": 0.05, "unit": "uM"}]


def test_metadata_selector_matches_same_plate_and_cell_line():
    rows = []
    for line in ("line-a", "line-b"):
        for drug in ("drug-a", "drug-b"):
            for replicate in range(2):
                rows.append({
                    "sample": f"{line}-{drug}-{replicate}", "plate": "plate1",
                    "cell_line_id": line, "drug": drug,
                    "drugname_drugconc": f"[('{drug}', 0.1, 'uM')]",
                    "n_cells": 20 + replicate,
                })
        for replicate in range(2):
            rows.append({
                "sample": f"{line}-control-{replicate}", "plate": "plate1",
                "cell_line_id": line, "drug": "DMSO_TF",
                "drugname_drugconc": "[]", "n_cells": 25,
            })
    frame = pd.DataFrame(rows)
    manifest = select_metadata_first_pilot(frame, conditions=3, min_cells=10)
    assert len({x["cell_line_id"] for x in manifest["conditions"]}) == 2
    assert all(x["treatment_sample_n"] == 2 and x["control_sample_n"] == 2 for x in manifest["conditions"])
    assert all(x["matched_plates"] == ["plate1"] for x in manifest["conditions"])


def test_control_without_plate_or_cell_line_match_is_not_selected():
    frame = pd.DataFrame([
        {"sample": "t1", "plate": "plate1", "cell_line_id": "line-a", "drug": "d1",
         "drugname_drugconc": "[]", "n_cells": 20},
        {"sample": "t2", "plate": "plate1", "cell_line_id": "line-a", "drug": "d1",
         "drugname_drugconc": "[]", "n_cells": 20},
        {"sample": "c1", "plate": "plate2", "cell_line_id": "line-a", "drug": "DMSO_TF",
         "drugname_drugconc": "[]", "n_cells": 20},
        {"sample": "c2", "plate": "plate2", "cell_line_id": "line-a", "drug": "DMSO_TF",
         "drugname_drugconc": "[]", "n_cells": 20},
    ])
    with pytest.raises(ValueError, match="cell lines"):
        select_metadata_first_pilot(frame, conditions=3, min_cells=10)


def test_partition_leakage_guard():
    clean = pd.DataFrame({"sample": ["a", "b"], "partition": ["train", "test"]})
    validate_group_partition(clean, group_column="sample")
    leaking = pd.DataFrame({"replicate": ["r1", "r1"], "partition": ["train", "test"]})
    with pytest.raises(ValueError, match="leakage"):
        validate_group_partition(leaking, group_column="replicate")
