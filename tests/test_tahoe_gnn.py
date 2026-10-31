import pytest


def test_gene_graph_model_cpu_forward_backward_and_checkpoint_reproducibility():
    torch = pytest.importorskip("torch")
    from scpolaseq_utrnet.tahoe_gnn import GeneGraphDeltaPredictor

    torch.manual_seed(7)
    model = GeneGraphDeltaPredictor(n_genes=5, context_dim=3, latent_dim=8, layers=2)
    basal = torch.rand(4, 5, requires_grad=True)
    context = torch.rand(4, 3)
    edges = torch.tensor([[0, 1, 2, 3, 4, 1], [1, 2, 3, 4, 0, 0]], dtype=torch.long)
    prediction = model(basal, context, edges)
    assert prediction.shape == (4, 5)
    prediction.square().mean().backward()
    assert basal.grad is not None and torch.isfinite(basal.grad).all()
    assert all(parameter.grad is not None for parameter in model.parameters())

    state = model.state_dict()
    clone = GeneGraphDeltaPredictor(n_genes=5, context_dim=3, latent_dim=8, layers=2)
    clone.load_state_dict(state)
    with torch.no_grad():
        assert torch.allclose(model(basal.detach(), context, edges), clone(basal.detach(), context, edges))
