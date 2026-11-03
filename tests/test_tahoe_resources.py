from scipy import sparse

from scpolaseq_utrnet.tahoe_resources import estimate_gnn_activation_bytes, sparse_resource_report


def test_sparse_resource_report_separates_payload_from_observed_peak():
    matrix = sparse.csr_matrix([[1, 0, 2], [0, 0, 3]])
    report = sparse_resource_report(matrix, elapsed_seconds=2.0, observed_peak_ram_bytes=1000)
    assert report["nnz"] == 3
    assert report["csr_payload_bytes"] == matrix.data.nbytes + matrix.indices.nbytes + matrix.indptr.nbytes
    assert report["csr_payload_is_ram_peak"] is False
    assert report["cells_per_second"] == 1.0
    assert report["observed_peak_ram_bytes"] == 1000


def test_gnn_activation_memory_estimate():
    result = estimate_gnn_activation_bytes(batch_size=2, n_genes=500, latent_dim=64, layers=2)
    assert result["activation_bytes_estimate"] == 2 * 500 * 64 * 4 * 4
