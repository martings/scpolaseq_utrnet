"""Generate a tiny self-contained UTRNet integration fixture."""
from __future__ import annotations

from pathlib import Path
import sys
import anndata as ad
import numpy as np
import pandas as pd
import pysam
import scipy.sparse as sp
import yaml


def make_fixture(root: Path) -> Path:
    upstream = root / "upstream"
    downstream = root / "downstream"
    for rel in ["reference", "labels", "pas_reference", "sierra_quant/cluster", "sierra_quant/cell_type"]:
        (upstream / rel).mkdir(parents=True, exist_ok=True)
    (downstream / "cell_labels").mkdir(parents=True, exist_ok=True)
    sequence = list("A" * 400)
    sequence[129:134] = list("ATTTA")  # retained only in distal plus-strand isoform
    sequence[209:214] = list("TAAAT")  # reverse-complement AUUUA on minus strand
    fasta = upstream / "reference/reference.genome.fa"
    fasta.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(fasta))
    gtf = upstream / "reference/reference.annotation.gtf"
    records = []
    def feature(kind, start, end, gene, transcript, strand):
        records.append(f'chr1\tfixture\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}.1"; transcript_id "{transcript}.1";\n')
    for kind, start, end in [("exon", 101, 140), ("CDS", 101, 110), ("UTR", 111, 140)]:
        feature(kind, start, end, "ENSGPLUS", "ENSTPLUS", "+")
    for kind, start, end in [("exon", 201, 240), ("CDS", 220, 240), ("UTR", 201, 219)]:
        feature(kind, start, end, "ENSGMINUS", "ENSTMINUS", "-")
    feature("exon", 120, 150, "ENSGOVERLAP", "ENSTOVERLAP", "-")
    gtf.write_text("".join(records))
    cells = pd.DataFrame([
        {"sample_id": "sample", "library_id": "sample", "barcode_raw": barcode,
         "barcode_corrected": barcode, "cell_id": f"sample:sample:{barcode}",
         "cluster_id": cluster, "cell_type": "unlabeled", "condition": "", "batch": "", "label_source": "fixture"}
        for barcode, cluster in [("BC1-1", "cluster_1"), ("BC2-1", "cluster_2")]])
    cells.to_csv(upstream / "labels/cell_annotations.tsv", sep="\t", index=False)
    obs = pd.DataFrame({"cell_id": cells.cell_id.to_list(), "celltype_corrected": ["A", "B"],
                        "qc_pass": [True, True], "cluster_id": ["cluster_1", "cluster_2"]}, index=cells.cell_id)
    var = pd.DataFrame({"gene_id": ["ENSGPLUS.1", "ENSGMINUS.1", "ENSGELAVL1.1", "ENSGZFP36.1"],
                        "gene_id_base": ["ENSGPLUS", "ENSGMINUS", "ENSGELAVL1", "ENSGZFP36"],
                        "gene_name": ["PLUS", "MINUS", "ELAVL1", "ZFP36"]}, index=["PLUS", "MINUS", "ELAVL1", "ZFP36"])
    counts = sp.csr_matrix(np.array([[5, 2, 2, 0], [8, 3, 4, 1]], dtype=np.int64))
    a = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    a.layers["counts"] = counts
    a.write_h5ad(downstream / "cell_labels/adata_celltyped.h5ad")
    pas = pd.DataFrame([
        ["pas:plus:125", "plus:125", "ENSGPLUS", "chr1", 125, 125, "+", "fixture"],
        ["pas:plus:140", "plus:140", "ENSGPLUS", "chr1", 140, 140, "+", "fixture"],
        ["pas:minus:215", "minus:215", "ENSGMINUS", "chr1", 215, 215, "-", "fixture"],
        ["pas:minus:205", "minus:205", "ENSGMINUS", "chr1", 205, 205, "-", "fixture"],
    ], columns=["pas_reference_id", "site_id", "gene_id", "chrom", "start", "end", "strand", "reference_source"])
    pas.to_csv(upstream / "pas_reference/pas_reference.tsv", sep="\t", index=False)
    cols = ["library_id", "group_level", "group_id", "site_id", "cell_barcode", "umi_count"]
    first = pd.DataFrame([["sample", "cluster", "cluster_1", "plus:125", "BC1-1", 8],
                  ["sample", "cluster", "cluster_1", "plus:140", "BC1-1", 2],
                  ["sample", "cluster", "cluster_1", "minus:215", "BC1-1", 3]], columns=cols)
    first.to_csv(upstream / "sierra_quant/cluster/sample.cluster.cluster_1.sierra_quant.tsv", sep="\t", index=False)
    second = pd.DataFrame([["sample", "cluster", "cluster_2", "plus:125", "BC2-1", 2],
                  ["sample", "cluster", "cluster_2", "plus:140", "BC2-1", 8],
                  ["sample", "cluster", "cluster_2", "minus:205", "BC2-1", 4]], columns=cols)
    second.to_csv(upstream / "sierra_quant/cluster/sample.cluster.cluster_2.sierra_quant.tsv", sep="\t", index=False)
    redundant = pd.concat([first, second], ignore_index=True)
    redundant["group_level"] = "cell_type"
    redundant["group_id"] = "unlabeled"
    redundant.to_csv(upstream / "sierra_quant/cell_type/sample.cell_type.unlabeled.sierra_quant.tsv", sep="\t", index=False)
    motif = root / "motifs.tsv"
    motif.write_text("motif_id\tregulator\tregulator_class\tpattern\tevidence_tier\tsource\tversion\turl\n"
                     "ARE\tELAVL1\tARE-BP\tAUUUA\tmotif_heuristic\tfixture\t1\t\n")
    mirna = root / "mirna.tsv"
    mirna.write_text("gene_id\tregulator\tchrom\tstart\tend\tstrand\tevidence_tier\tsource\tversion\turl\n"
                     "ENSGPLUS\thsa-miR-fixture\tchr1\t130\t134\t+\tcomputational_prediction\tfixture\t1\t\n")
    config = {"contract": "utrnet-input/v0.1", "sample_id": "sample", "gene": "PLUS",
              "upstream": str(upstream), "downstream": str(downstream), "output": str(root / "outputs"),
              "min_population_cells": 1, "resources": {"motifs": str(motif), "mirna": str(mirna)}}
    target = root / "config.yaml"
    target.write_text(yaml.safe_dump(config))
    return target


if __name__ == "__main__":
    print(make_fixture(Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/utrnet-fixture").resolve()))
