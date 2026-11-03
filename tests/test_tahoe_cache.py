import numpy as np
import pytest
from scipy import sparse

from scpolaseq_utrnet.tahoe_cache import write_sparse_cache


def test_budgeted_sparse_cache_round_trip(tmp_path):
    matrix = sparse.csr_matrix(np.array([[1, 0], [0, 2]], dtype=np.int64))
    path, sidecar = write_sparse_cache(
        matrix, ["cell-a", "cell-b"], ["ENSG1", "ENSG2"],
        tmp_path / "pilot.npz", max_bytes=100_000,
    )
    loaded = sparse.load_npz(path)
    record = __import__("json").loads(sidecar.read_text())
    assert (loaded != matrix).nnz == 0
    assert record["raw_count_sum"] == 3
    assert record["missing_modalities"] == ["PAS usage", "timepoint"]


def test_cache_rejects_budget_overflow_and_identifier_mismatch(tmp_path):
    matrix = sparse.csr_matrix(np.eye(10, dtype=np.int64))
    with pytest.raises(ValueError, match="exceeding"):
        write_sparse_cache(matrix, [str(i) for i in range(10)], [str(i) for i in range(10)],
                           tmp_path / "large.npz", max_bytes=1)
    with pytest.raises(ValueError, match="dimensions"):
        write_sparse_cache(matrix, ["cell"], ["gene"], tmp_path / "bad.npz", max_bytes=100_000)
