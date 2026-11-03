"""Small optional gene-graph delta predictor used by the Tahoe benchmark.

This dependency-light message-passing stub uses PyTorch only and does not
claim that graph edges are measured molecular interactions.
"""
from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class GeneGraphDeltaPredictor(nn.Module):
    """Predict one endpoint delta per gene from basal expression and context."""

    def __init__(self, n_genes: int, context_dim: int, latent_dim: int = 64, layers: int = 2):
        super().__init__()
        if n_genes < 1 or context_dim < 1 or latent_dim < 1 or layers < 1:
            raise ValueError("model dimensions and layer count must be positive")
        self.n_genes = n_genes
        self.gene_embedding = nn.Embedding(n_genes, latent_dim)
        self.expression_projection = nn.Linear(1, latent_dim)
        self.context_projection = nn.Linear(context_dim, latent_dim)
        self.self_layers = nn.ModuleList(nn.Linear(latent_dim, latent_dim) for _ in range(layers))
        self.neighbor_layers = nn.ModuleList(nn.Linear(latent_dim, latent_dim, bias=False) for _ in range(layers))
        self.output = nn.Linear(latent_dim, 1)

    def forward(self, basal_log_counts: torch.Tensor, context: torch.Tensor, edge_index: torch.Tensor) -> torch.Tensor:
        if basal_log_counts.ndim != 2 or basal_log_counts.shape[1] != self.n_genes:
            raise ValueError("basal_log_counts must have shape [batch, n_genes]")
        if context.ndim != 2 or context.shape[0] != basal_log_counts.shape[0]:
            raise ValueError("context must have shape [batch, context_features]")
        if edge_index.ndim != 2 or edge_index.shape[0] != 2:
            raise ValueError("edge_index must have shape [2, n_edges]")
        batch_size = basal_log_counts.shape[0]
        gene_ix = torch.arange(self.n_genes, device=basal_log_counts.device)
        nodes = self.gene_embedding(gene_ix).unsqueeze(0)
        nodes = nodes + self.expression_projection(basal_log_counts.unsqueeze(-1))
        nodes = nodes + self.context_projection(context).unsqueeze(1)
        src, dst = edge_index[0].long(), edge_index[1].long()
        if src.numel() and (int(torch.max(torch.stack([src.max(), dst.max()]))) >= self.n_genes or int(torch.min(torch.stack([src.min(), dst.min()]))) < 0):
            raise ValueError("edge index outside gene range")
        degree = torch.bincount(dst, minlength=self.n_genes).clamp_min(1).to(nodes.dtype).view(1, -1, 1)
        for self_layer, neighbor_layer in zip(self.self_layers, self.neighbor_layers):
            aggregate = torch.zeros_like(nodes)
            if src.numel():
                aggregate.index_add_(1, dst, nodes.index_select(1, src))
            nodes = F.gelu(self_layer(nodes) + neighbor_layer(aggregate / degree))
        return self.output(nodes).squeeze(-1)
