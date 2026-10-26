"""Versioned motif and imported miRNA evidence, separate from measured occupancy."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from .contracts import read_tsv


SITE_COLUMNS = ["regulatory_site_id", "isoform_id", "site_id", "motif_id", "regulator",
                "regulator_class", "isoform_start", "isoform_end", "chrom", "genomic_start",
                "genomic_end", "strand", "evidence_tier", "source", "version", "source_url",
                "pubmed_id", "assay", "cell_context", "matrix_id"]


def _mature_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper().replace("T", "U")))
            name, chunks = line[1:].split()[0], []
        elif name is None:
            raise ValueError("miRNA FASTA has sequence before a header")
        else:
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper().replace("T", "U")))
    if not records or len({name for name, _ in records}) != len(records):
        raise ValueError("miRNA FASTA is empty or has duplicate IDs")
    for name, sequence in records:
        if len(sequence) < 8 or set(sequence) - set("ACGU"):
            raise ValueError(f"invalid mature miRNA sequence: {name}")
    return records


def scan_mirna_seeds(isoforms: pd.DataFrame, coord_maps: dict[str, list[int]],
                     sites: pd.DataFrame, resource: Path | None) -> pd.DataFrame:
    """Predict exact 7mer-m8 complements to mature miRNA bases 2–8.

    A seed match is a sequence prior; it is not a validated target or measured
    miRNA expression. Sites crossing a splice junction are skipped because the
    current output contract has only one contiguous genomic interval.
    """
    if resource is None:
        return pd.DataFrame(columns=SITE_COLUMNS)
    complement = str.maketrans("ACGU", "UGCA")
    seeds = [(name, sequence[1:8].translate(complement)[::-1])
             for name, sequence in _mature_fasta(resource)]
    site_lookup = sites.set_index("site_id")
    rows = []
    for iso in isoforms.itertuples(index=False):
        if not iso.sequence:
            continue
        sequence = iso.sequence.upper().replace("T", "U")
        coords = coord_maps[iso.isoform_id]
        pas = site_lookup.loc[iso.site_id]
        for name, target in seeds:
            for i in range(len(sequence) - 6):
                if sequence[i:i + 7] != target:
                    continue
                region = coords[i:i + 7]
                if len(region) != 7 or sorted(region) != list(range(min(region), max(region) + 1)):
                    continue
                rows.append({"regulatory_site_id": f"{iso.isoform_id}:mirna_seed:{name}:{i}",
                             "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                             "motif_id": f"mirna_seed:{name}:2-8", "regulator": name,
                             "regulator_class": "miRNA", "isoform_start": i,
                             "isoform_end": i + 7, "chrom": pas.chrom,
                             "genomic_start": min(region), "genomic_end": max(region),
                             "strand": pas.strand, "evidence_tier": "computational_prediction",
                             "source": "MirGeneDB mature sequence; 7mer-m8 seed scan",
                             "version": "3.0", "source_url": "https://www.mirgenedb.org/fasta/hsa?mat=1",
                             "pubmed_id": "", "assay": "in_silico_exact_seed_match",
                             "cell_context": "unmeasured_in_PBMC", "matrix_id": ""})
    return pd.DataFrame(rows, columns=SITE_COLUMNS)


def scan_motifs(isoforms: pd.DataFrame, coord_maps: dict[str, list[int]],
                sites: pd.DataFrame, resource: Path | None) -> pd.DataFrame:
    if resource is None:
        return pd.DataFrame(columns=SITE_COLUMNS)
    motifs = read_tsv(resource)
    required = {"motif_id", "regulator", "regulator_class", "pattern", "evidence_tier", "source", "version", "url"}
    if not required.issubset(motifs.columns):
        raise ValueError(f"motif resource lacks {sorted(required - set(motifs.columns))}")
    if motifs.motif_id.duplicated().any():
        raise ValueError("duplicate motif ID")
    site_lookup = sites.set_index("site_id")
    rows = []
    for iso in isoforms.itertuples(index=False):
        if not iso.sequence:
            continue
        coords = coord_maps[iso.isoform_id]
        pas = site_lookup.loc[iso.site_id]
        for motif in motifs.itertuples(index=False):
            pattern = motif.pattern.upper().replace("T", "U")
            if not pattern or set(pattern) - set("ACGU"):
                raise ValueError(f"unsupported motif alphabet: {motif.motif_id}")
            seq = iso.sequence.upper().replace("T", "U")
            for i in range(len(seq) - len(pattern) + 1):
                if seq[i:i + len(pattern)] != pattern:
                    continue
                region = coords[i:i + len(pattern)]
                rows.append({"regulatory_site_id": f"{iso.isoform_id}:{motif.motif_id}:{i}",
                             "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                             "motif_id": motif.motif_id, "regulator": motif.regulator,
                             "regulator_class": motif.regulator_class,
                             "isoform_start": i, "isoform_end": i + len(pattern),
                             "chrom": pas.chrom, "genomic_start": min(region),
                             "genomic_end": max(region), "strand": pas.strand,
                             "evidence_tier": motif.evidence_tier, "source": motif.source,
                             "version": motif.version, "source_url": motif.url,
                             "pubmed_id": getattr(motif, "pubmed_id", ""),
                             "assay": getattr(motif, "assay", ""),
                             "cell_context": getattr(motif, "cell_context", ""),
                             "matrix_id": getattr(motif, "matrix_id", "")})
    return pd.DataFrame(rows, columns=SITE_COLUMNS)


def import_mirna_sites(isoforms: pd.DataFrame, coord_maps: dict[str, list[int]],
                       resource: Path | None) -> pd.DataFrame:
    if resource is None:
        return pd.DataFrame(columns=SITE_COLUMNS)
    frame = read_tsv(resource)
    required = {"gene_id", "regulator", "chrom", "start", "end", "strand", "evidence_tier", "source", "version", "url"}
    if not required.issubset(frame.columns):
        raise ValueError(f"miRNA resource lacks {sorted(required - set(frame.columns))}")
    rows = []
    for iso in isoforms.itertuples(index=False):
        if not iso.sequence:
            continue
        coords = coord_maps[iso.isoform_id]
        for hit in frame.itertuples(index=False):
            if hit.gene_id.split(".")[0] != iso.gene_id:
                continue
            start, end = int(hit.start), int(hit.end)
            matches = [i for i, c in enumerate(coords) if start <= c <= end]
            if len(matches) != end - start + 1 or matches != list(range(min(matches), max(matches) + 1)):
                continue
            i = min(matches)
            rows.append({"regulatory_site_id": f"{iso.isoform_id}:mirna:{hit.regulator}:{i}",
                         "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                         "motif_id": f"mirna:{hit.regulator}:{hit.chrom}:{start}-{end}",
                         "regulator": hit.regulator, "regulator_class": "miRNA",
                         "isoform_start": i, "isoform_end": max(matches) + 1,
                         "chrom": hit.chrom, "genomic_start": start, "genomic_end": end,
                         "strand": hit.strand, "evidence_tier": hit.evidence_tier,
                         "source": hit.source, "version": hit.version, "source_url": hit.url})
    return pd.DataFrame(rows, columns=SITE_COLUMNS)
