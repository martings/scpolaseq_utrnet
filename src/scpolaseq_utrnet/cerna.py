"""Transcriptome-scale competing endogenous RNA (ceRNA) and RBP-sponge network modeling."""
from __future__ import annotations

from dataclasses import replace
import json
import logging
from pathlib import Path
from typing import Any
import networkx as nx
import numpy as np
import pandas as pd

from .config import Config
from .ingest import ingest_gene
from .isoforms import reconstruct_isoforms
from .motifs import scan_motifs

logger = logging.getLogger(__name__)


def get_top_prioritized_genes(config: Config, top_n: int = 15) -> list[str]:
    """Retrieve the top prioritized APA genes from downstream priority table."""
    if not config.apa_priority_table.exists():
        logger.warning("apa_priority_table not found at %s", config.apa_priority_table)
        return [config.gene]

    df = pd.read_csv(config.apa_priority_table, sep="\t")
    if "is_robust" in df.columns:
        df = df[df.is_robust.astype(bool)]
    if "basis_span_bp" in df.columns:
        df = df[df.basis_span_bp <= config.max_utr_span]
    if df.empty:
        return [config.gene]

    # Rank by maximum priority score per gene
    ranked = (
        df.groupby("gene_name")["priority_score"]
        .max()
        .sort_values(ascending=False)
        .index.tolist()
    )
    # Clean gene symbols (handle composite IDs like "GENEA|GENEB")
    genes = []
    for g in ranked:
        clean = str(g).split("|")[0].strip()
        if clean and clean not in genes:
            genes.append(clean)
        if len(genes) >= top_n:
            break
    if config.gene not in genes and config.gene:
        genes.insert(0, config.gene)
    return genes[:top_n]


def extract_gene_sites(
    config: Config, gene: str, cache_dir: Path | None = None
) -> pd.DataFrame:
    """Extract binding sites per PAS and regulator for a single gene."""
    # Check cache directory if regulatory_sites.csv already exists
    if cache_dir is not None:
        cand_files = list(cache_dir.glob(f"*/regulatory_sites.csv"))
        # We can also check if subfolder matches gene or gene_id
        for cf in cand_files:
            try:
                sub = pd.read_csv(cf)
                if not sub.empty and "regulator" in sub.columns and "site_id" in sub.columns:
                    # Check if gene matches
                    if "isoform_id" in sub.columns and gene in str(sub.isoform_id.iloc[0]):
                        return sub[["site_id", "regulator", "regulator_class", "motif_id"]]
            except Exception:
                pass

    try:
        gene_cfg = replace(config, gene=gene)
        tables = ingest_gene(gene_cfg)
        if tables["pas_sites"].empty:
            return pd.DataFrame(columns=["site_id", "regulator", "regulator_class", "motif_id"])

        isoforms, maps = reconstruct_isoforms(
            tables["pas_sites"],
            gene_cfg.reference_gtf,
            gene_cfg.reference_fasta,
            max_span=gene_cfg.max_utr_span,
        )
        if isoforms.empty or (isoforms.sequence_length <= 0).all():
            logger.warning("Gene %s lacks valid compatible 3'UTR models; skipping from ceRNA cohort", gene)
            return pd.DataFrame(columns=["site_id", "regulator", "regulator_class", "motif_id"])

        sites = scan_motifs(isoforms, maps, tables["pas_sites"], gene_cfg.motif_resource)
        if sites.empty:
            return pd.DataFrame(columns=["site_id", "regulator", "regulator_class", "motif_id"])

        return sites[["site_id", "regulator", "regulator_class", "motif_id"]]
    except Exception as exc:
        logger.warning("Failed to extract sites for gene %s: %s", gene, exc)
        return pd.DataFrame(columns=["site_id", "regulator", "regulator_class", "motif_id"])


def compute_cerna_load(
    config: Config,
    genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute transcriptome-scale binding site load and effective sites per gene."""
    gene_site_records = []
    all_profiles = []

    for g in genes:
        try:
            gene_cfg = replace(config, gene=g)
            tables = ingest_gene(gene_cfg)
            sites = extract_gene_sites(config, g, cache_dir=config.output)
            if sites.empty or tables["apa_profiles"].empty:
                continue

            # Count sites per site_id and regulator
            site_counts = (
                sites.groupby(["site_id", "regulator", "regulator_class"])
                .size()
                .reset_index(name="n_sites")
            )
            site_counts["gene"] = g

            gene_site_records.append(site_counts)

            prof = tables["apa_profiles"].copy()
            prof["gene"] = g
            all_profiles.append(prof)
        except Exception as exc:
            logger.warning("Skipping gene %s in ceRNA calculation: %s", g, exc)

    if not gene_site_records or not all_profiles:
        empty_sites = pd.DataFrame(columns=["population", "gene", "regulator", "effective_sites", "gene_expression", "gene_load"])
        empty_load = pd.DataFrame(columns=["population", "regulator", "total_load"])
        empty_counts = pd.DataFrame(columns=["gene", "site_id", "regulator", "regulator_class", "n_sites"])
        return empty_sites, empty_load, empty_counts

    combined_sites = pd.concat(gene_site_records, ignore_index=True)
    combined_profiles = pd.concat(all_profiles, ignore_index=True)

    # For each population, gene, and regulator, compute:
    # Effective sites N_r(g, c) = sum_p (Usage_p(c) * Sites_p(r))
    merged = combined_profiles.merge(
        combined_sites, on=["gene", "site_id"], how="inner"
    )
    merged["weighted_sites"] = merged["pas_proportion"].fillna(0.0) * merged["n_sites"]
    if "total_gene_umi" in merged.columns and "n_cells" in merged.columns:
        merged["gene_expression"] = merged["total_gene_umi"] / merged["n_cells"].clip(lower=1)
    elif "mean_gene_umi" in merged.columns:
        merged["gene_expression"] = merged["mean_gene_umi"]
    else:
        merged["gene_expression"] = 1.0

    effective_df = (
        merged.groupby(["celltype_corrected", "gene", "regulator", "regulator_class"], as_index=False)
        .agg(
            effective_sites=("weighted_sites", "sum"),
            gene_expression=("gene_expression", "first"),
        )
        .rename(columns={"celltype_corrected": "population"})
    )
    effective_df["gene_load"] = effective_df["effective_sites"] * effective_df["gene_expression"]

    # Total cellular load for each regulator in each population
    load_df = (
        effective_df.groupby(["population", "regulator", "regulator_class"], as_index=False)["gene_load"]
        .sum()
        .rename(columns={"gene_load": "total_load"})
    )

    return effective_df, load_df, combined_sites


def compute_cerna_coupling(
    effective_df: pd.DataFrame,
    load_df: pd.DataFrame,
    population: str,
    min_coupling: float = 0.02,
) -> pd.DataFrame:
    """Compute pairwise ceRNA cross-talk coupling scores in a specific population."""
    eff = effective_df[effective_df.population == population].copy()
    loads = load_df[load_df.population == population].set_index("regulator")["total_load"].to_dict()

    if eff.empty:
        return pd.DataFrame(columns=[
            "population", "gene_a", "gene_b", "coupling_score",
            "shared_regulators_count", "shared_regulators_list"
        ])

    genes = sorted(eff.gene.unique())
    # Pre-map gene loads: gene -> dict of regulator -> load
    gene_reg_loads: dict[str, dict[str, float]] = {g: {} for g in genes}
    for r in eff.itertuples(index=False):
        gene_reg_loads[r.gene][r.regulator] = float(r.gene_load)

    # Compute self-coupling for normalization
    self_coupling: dict[str, float] = {}
    for g in genes:
        val = sum(
            (load_val ** 2) / (loads.get(reg, 1.0) + 1e-6)
            for reg, load_val in gene_reg_loads[g].items()
        )
        self_coupling[g] = max(val, 1e-9)

    records = []
    for i, g1 in enumerate(genes):
        for g2 in genes[i + 1:]:
            common_regs = set(gene_reg_loads[g1].keys()).intersection(gene_reg_loads[g2].keys())
            if not common_regs:
                continue

            raw_coupling = sum(
                (gene_reg_loads[g1][r] * gene_reg_loads[g2][r]) / (loads.get(r, 1.0) + 1e-6)
                for r in common_regs
            )
            norm_coupling = raw_coupling / np.sqrt(self_coupling[g1] * self_coupling[g2])

            if norm_coupling >= min_coupling:
                records.append({
                    "population": population,
                    "gene_a": g1,
                    "gene_b": g2,
                    "coupling_score": float(norm_coupling),
                    "shared_regulators_count": len(common_regs),
                    "shared_regulators_list": ";".join(sorted(common_regs)),
                })

    return pd.DataFrame(records)


def build_cerna_network(
    coupling_df: pd.DataFrame,
    effective_df: pd.DataFrame,
    load_df: pd.DataFrame,
    population: str,
) -> nx.Graph:
    """Construct ceRNA trans-competition network graph for a given cell population."""
    graph = nx.Graph(population=population, model="utrnet-v0.2-cerna")

    eff_sub = effective_df[effective_df.population == population]
    load_sub = load_df[load_df.population == population]

    # Add gene nodes
    for gene, grp in eff_sub.groupby("gene"):
        expr = float(grp["gene_expression"].iloc[0])
        total_sites = float(grp["effective_sites"].sum())
        total_load = float(grp["gene_load"].sum())
        graph.add_node(
            f"gene:{gene}",
            kind="gene",
            label=gene,
            expression=expr,
            effective_sites=total_sites,
            cellular_site_load=total_load,
        )

    # Add regulator nodes
    for r in load_sub.itertuples(index=False):
        graph.add_node(
            f"regulator:{r.regulator}",
            kind=r.regulator_class,
            label=r.regulator,
            total_cellular_load=float(r.total_load),
        )

    # Add gene -- regulator titration edges (top 50 strongest titration edges)
    top_titrations = eff_sub.sort_values("gene_load", ascending=False).head(100)
    for r in top_titrations.itertuples(index=False):
        graph.add_edge(
            f"gene:{r.gene}",
            f"regulator:{r.regulator}",
            kind="titrates_regulator",
            weight=float(r.gene_load),
            effective_sites=float(r.effective_sites),
        )

    # Add pairwise ceRNA cross-talk edges
    for r in coupling_df.itertuples(index=False):
        graph.add_edge(
            f"gene:{r.gene_a}",
            f"gene:{r.gene_b}",
            kind="cerna_crosstalk",
            weight=float(r.coupling_score),
            shared_regulators=r.shared_regulators_list,
            n_shared=int(r.shared_regulators_count),
        )

    return graph


def export_cerna_analysis(
    config: Config,
    genes: list[str] | None = None,
    top_n: int = 15,
    output_dir: Path | None = None,
) -> dict[str, Path]:
    """Execute end-to-end ceRNA and RBP-sponge analysis across prioritized APA targets."""
    out = output_dir or (config.output / "cerna")
    out.mkdir(parents=True, exist_ok=True)

    target_genes = genes or get_top_prioritized_genes(config, top_n=top_n)
    logger.info("Executing ceRNA network analysis on %d genes: %s", len(target_genes), target_genes)

    effective_df, load_df, site_counts = compute_cerna_load(config, target_genes)

    effective_path = out / "cerna_gene_effective_sites.csv"
    load_path = out / "cerna_load_by_celltype.csv"
    effective_df.to_csv(effective_path, index=False)
    load_df.to_csv(load_path, index=False)

    all_couplings = []
    populations = sorted(effective_df.population.unique()) if not effective_df.empty else []

    graph_dir = out / "graphs"
    graph_dir.mkdir(parents=True, exist_ok=True)

    for pop in populations:
        coup = compute_cerna_coupling(effective_df, load_df, pop)
        if not coup.empty:
            all_couplings.append(coup)
        graph = build_cerna_network(coup, effective_df, load_df, pop)
        slug = "".join(c.lower() if c.isalnum() else "_" for c in pop).strip("_")
        (graph_dir / f"{slug}_cerna.cytoscape.json").write_text(
            json.dumps(nx.cytoscape_data(graph), indent=2)
        )
        nx.write_graphml(graph, graph_dir / f"{slug}_cerna.graphml")

    coupling_matrix_path = out / "cerna_coupling_matrix.csv"
    combined_coupling = pd.concat(all_couplings, ignore_index=True) if all_couplings else pd.DataFrame(
        columns=["population", "gene_a", "gene_b", "coupling_score", "shared_regulators_count", "shared_regulators_list"]
    )
    combined_coupling.to_csv(coupling_matrix_path, index=False)

    return {
        "output_dir": out,
        "effective_sites": effective_path,
        "cellular_load": load_path,
        "coupling_matrix": coupling_matrix_path,
        "graphs": graph_dir,
    }
