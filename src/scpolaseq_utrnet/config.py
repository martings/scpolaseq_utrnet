"""Versioned configuration and path resolution."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import os
import yaml


@dataclass(frozen=True)
class Config:
    upstream: Path
    downstream: Path
    output: Path
    sample_id: str
    gene: str
    min_population_cells: int = 15
    min_pas_umi: int = 5
    max_utr_span: int = 50000
    fold_window: int = 80
    fold_span: int = 40
    fold_unpaired: int = 5
    motif_resource: Path | None = None
    mirna_resource: Path | None = None
    mirna_seed_resource: Path | None = None

    @property
    def curated_h5ad(self) -> Path:
        return self.downstream / "cell_labels/adata_celltyped.h5ad"

    @property
    def apa_priority_table(self) -> Path:
        return self.downstream / "apa/celltype_specific_APA_priority_table.tsv"

    @property
    def si_magnetization_table(self) -> Path:
        return self.downstream / "apa/si_and_magnetization_by_celltype.tsv"

    @property
    def reference_fasta(self) -> Path:
        return self.upstream / "reference/reference.genome.fa"

    @property
    def reference_gtf(self) -> Path:
        return self.upstream / "reference/reference.annotation.gtf"


def load_config(path: str | Path) -> Config:
    p = Path(path).expanduser().resolve()
    data = yaml.safe_load(p.read_text())
    if data.get("contract") not in {"utrnet-input/v0.1", "utrnet-input/v0.2"}:
        raise ValueError("config contract must be utrnet-input/v0.1 or utrnet-input/v0.2")

    def locate(key: str, fallback: str | None = None) -> Path:
        configured = str(data.get(key, fallback))
        raw = os.path.expandvars(configured)
        if "$" in raw:
            raise ValueError(f"unresolved environment variable in {key}: {configured}")
        candidate = Path(raw).expanduser()
        return (candidate if candidate.is_absolute() else p.parent / candidate).resolve()

    resources = data.get("resources", {})
    def resource(key: str) -> Path | None:
        raw = resources.get(key)
        if not raw:
            return None
        candidate = Path(os.path.expandvars(str(raw))).expanduser()
        return (candidate if candidate.is_absolute() else p.parent / candidate).resolve()

    return Config(
        upstream=locate("upstream"), downstream=locate("downstream"),
        output=locate("output", "../outputs/pbmc1k"),
        sample_id=str(data["sample_id"]), gene=str(data.get("gene", "SAT1")),
        min_population_cells=int(data.get("min_population_cells", 15)),
        min_pas_umi=int(data.get("min_pas_umi", 5)),
        max_utr_span=int(data.get("max_utr_span", 50000)),
        fold_window=int(data.get("fold_window", 80)),
        fold_span=int(data.get("fold_span", 40)),
        fold_unpaired=int(data.get("fold_unpaired", 5)),
        motif_resource=resource("motifs"), mirna_resource=resource("mirna"),
        mirna_seed_resource=resource("mirna_seeds"),
    )
