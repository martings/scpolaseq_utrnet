"""End-to-end frozen-output graph build."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pandas as pd

from .config import Config
from .ingest import ingest_gene, write_ingestion
from .isoforms import reconstruct_isoforms
from .motifs import scan_motifs, import_mirna_sites, scan_mirna_seeds
from .folding import annotate_accessibility, common_region_accessibility
from .regulators import regulator_abundance
from .graphs import score_interactions, compare_populations, export_graphs, site_retention, population_metrics


def _fingerprint(path: Path) -> dict:
    stat = path.stat()
    digest = None
    if stat.st_size <= 100_000_000:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                h.update(block)
        digest = h.hexdigest()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns, "sha256": digest}


def build_gene_graphs(config: Config) -> dict[str, Path]:
    tables = ingest_gene(config)
    gene_id = str(tables["pas_sites"].gene_id.iloc[0]).split(".")[0]
    directory = config.output / gene_id
    write_ingestion(tables, directory)
    isoforms, maps = reconstruct_isoforms(tables["pas_sites"], config.reference_gtf,
                                          config.reference_fasta, config.max_utr_span)
    isoforms.to_csv(directory / "utr_isoforms.csv", index=False)
    motif_sites = scan_motifs(isoforms, maps, tables["pas_sites"], config.motif_resource)
    mirna_sites = import_mirna_sites(isoforms, maps, config.mirna_resource)
    seed_sites = scan_mirna_seeds(isoforms, maps, tables["pas_sites"], config.mirna_seed_resource)
    regulatory_sites = pd.concat([motif_sites, mirna_sites, seed_sites], ignore_index=True)
    regulatory_sites = annotate_accessibility(regulatory_sites, isoforms, directory / "fold_cache",
        window=config.fold_window, span=config.fold_span, unpaired=config.fold_unpaired)
    regulatory_sites.to_csv(directory / "regulatory_sites.csv", index=False)
    site_retention(isoforms, regulatory_sites).to_csv(directory / "site_retention.csv", index=False)
    common_region_accessibility(regulatory_sites).to_csv(directory / "common_region_accessibility.csv", index=False)
    evidence = regulatory_sites[["motif_id", "regulator", "evidence_tier", "source", "version", "source_url",
                                 "pubmed_id", "assay", "cell_context", "matrix_id"]].drop_duplicates().copy()
    evidence = evidence.rename(columns={"motif_id": "accession", "source": "database", "source_url": "url"})
    evidence["site_detection_method"] = "sequence_scan_or_import"
    evidence.to_csv(directory / "evidence.csv", index=False)
    abundance = regulator_abundance(config, regulatory_sites, tables["cells"])
    abundance.to_csv(directory / "regulator_abundance.csv", index=False)
    scores = score_interactions(tables["apa_profiles"], isoforms, regulatory_sites,
                                abundance, config.min_population_cells)
    scores.to_csv(directory / "interaction_scores.csv", index=False)
    population_metrics(scores).to_csv(directory / "population_metrics.csv", index=False)
    rewiring = compare_populations(scores)
    rewiring.to_csv(directory / "rewiring.csv", index=False)
    export_graphs(scores, tables["apa_profiles"], isoforms, regulatory_sites, directory / "graphs")
    sources = [config.curated_h5ad, config.upstream / "pas_reference/pas_reference.tsv",
               config.reference_gtf, config.reference_fasta]
    if config.motif_resource:
        sources.append(config.motif_resource)
    if config.mirna_resource:
        sources.append(config.mirna_resource)
    if config.mirna_seed_resource:
        sources.append(config.mirna_seed_resource)
    metadata = {"contract": "utrnet-input/v0.1", "model": "utrnet-v0.1",
                "sample_id": config.sample_id, "gene_id": gene_id,
                "source_execution": str(config.upstream), "sources": [_fingerprint(p) for p in sources],
                "folding_parameters": {"window": config.fold_window, "span": config.fold_span,
                                       "unpaired": config.fold_unpaired},
                "limitations": ["single biological sampling context", "interaction scores are predictions, not occupancy",
                                "miRNA abundance is unmeasured without small-RNA input",
                                "unresolved or ambiguous transcript/PAS associations remain flagged"]}
    (directory / "provenance.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return {"directory": directory, "scores": directory / "interaction_scores.csv",
            "graphs": directory / "graphs"}
