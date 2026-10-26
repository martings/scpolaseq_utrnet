"""Versioned motif and imported miRNA evidence, separate from measured occupancy."""
from __future__ import annotations

from pathlib import Path
import re
import pandas as pd

from .contracts import read_tsv


SITE_COLUMNS = ["regulatory_site_id", "isoform_id", "site_id", "motif_id", "regulator",
                "regulator_class", "isoform_start", "isoform_end", "chrom", "genomic_start",
                "genomic_end", "strand", "evidence_tier", "source", "version", "source_url",
                "pubmed_id", "assay", "cell_context", "matrix_id"]

IUPAC_RNA_ALPHABET = set("ACGU" + "RYSWKMBDHVN")

IUPAC_TO_REGEX = {
    "A": "A", "C": "C", "G": "G", "U": "U",
    "R": "[AG]", "Y": "[CU]", "S": "[GC]", "W": "[AU]",
    "K": "[GU]", "M": "[AC]", "B": "[CGU]", "D": "[AGU]",
    "H": "[ACU]", "V": "[ACG]", "N": "[ACGU]"
}

SEED_WEIGHTS = {
    "8mer": 1.0,
    "7mer-m8": 0.8,
    "7mer-A1": 0.6,
    "6mer": 0.3,
}


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
                     sites: pd.DataFrame, resource: Path | None,
                     seed_types: tuple[str, ...] = ("7mer-m8",)) -> pd.DataFrame:
    """Predict calibrated canonical complements to mature miRNA seeds.

    Supported seed types:
      - 8mer: bases 2–8 match + A opposite position 1 (weight 1.0)
      - 7mer-m8: exact complement to bases 2–8 (weight 0.8)
      - 7mer-A1: bases 2–7 match + A opposite position 1 (weight 0.6)
      - 6mer: exact complement to bases 2–7 (weight 0.3)

    A seed match is a sequence prior; it is not a validated target or measured
    miRNA expression. Sites crossing a splice junction are skipped because the
    current output contract has only one contiguous genomic interval.
    """
    if resource is None:
        return pd.DataFrame(columns=SITE_COLUMNS)
    complement = str.maketrans("ACGU", "UGCA")
    mature_entries = _mature_fasta(resource)
    site_lookup = sites.set_index("site_id")
    rows = []

    for iso in isoforms.itertuples(index=False):
        if not iso.sequence:
            continue
        sequence = iso.sequence.upper().replace("T", "U")
        coords = coord_maps[iso.isoform_id]
        pas = site_lookup.loc[iso.site_id]

        for name, mat_seq in mature_entries:
            # 7mer-m8 target: complement of bases 2-8 (1-indexed 2 to 8, 0-indexed 1 to 7)
            target_7m8 = mat_seq[1:8].translate(complement)[::-1]
            # 6mer target: complement of bases 2-7
            target_6mer = mat_seq[1:7].translate(complement)[::-1]

            if "7mer-m8" in seed_types:
                t_len = len(target_7m8)
                for i in range(len(sequence) - t_len + 1):
                    if sequence[i:i + t_len] != target_7m8:
                        continue
                    region = coords[i:i + t_len]
                    if len(region) != t_len or sorted(region) != list(range(min(region), max(region) + 1)):
                        continue
                    rows.append({"regulatory_site_id": f"{iso.isoform_id}:mirna_seed:{name}:{i}",
                                 "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                                 "motif_id": f"mirna_seed:{name}:2-8", "regulator": name,
                                 "regulator_class": "miRNA", "isoform_start": i,
                                 "isoform_end": i + t_len, "chrom": pas.chrom,
                                 "genomic_start": min(region), "genomic_end": max(region),
                                 "strand": pas.strand, "evidence_tier": "computational_prediction",
                                 "source": "MirGeneDB mature sequence; 7mer-m8 seed scan",
                                 "version": "3.0", "source_url": "https://www.mirgenedb.org/fasta/hsa?mat=1",
                                 "pubmed_id": "26432245", "assay": "in_silico_canonical_7mer_m8",
                                 "cell_context": "unmeasured_in_PBMC", "matrix_id": "7mer-m8"})

            if "8mer" in seed_types:
                # 8mer: target_7m8 followed by 'A' at 3' end in mRNA (opposite miRNA base 1)
                target_8mer = target_7m8 + "A"
                t_len = len(target_8mer)
                for i in range(len(sequence) - t_len + 1):
                    if sequence[i:i + t_len] != target_8mer:
                        continue
                    region = coords[i:i + t_len]
                    if len(region) != t_len or sorted(region) != list(range(min(region), max(region) + 1)):
                        continue
                    rows.append({"regulatory_site_id": f"{iso.isoform_id}:mirna_seed_8mer:{name}:{i}",
                                 "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                                 "motif_id": f"mirna_seed:{name}:8mer", "regulator": name,
                                 "regulator_class": "miRNA", "isoform_start": i,
                                 "isoform_end": i + t_len, "chrom": pas.chrom,
                                 "genomic_start": min(region), "genomic_end": max(region),
                                 "strand": pas.strand, "evidence_tier": "computational_prediction",
                                 "source": "MirGeneDB mature sequence; 8mer seed scan",
                                 "version": "3.0", "source_url": "https://www.mirgenedb.org/fasta/hsa?mat=1",
                                 "pubmed_id": "26432245", "assay": "in_silico_canonical_8mer",
                                 "cell_context": "unmeasured_in_PBMC", "matrix_id": "8mer"})

    return pd.DataFrame(rows, columns=SITE_COLUMNS)


def scan_motifs(isoforms: pd.DataFrame, coord_maps: dict[str, list[int]],
                sites: pd.DataFrame, resource: Path | None) -> pd.DataFrame:
    """Scan 3'UTR isoforms for known RNA binding protein (RBP) motifs using IUPAC degeneracies."""
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
        seq = iso.sequence.upper().replace("T", "U")

        for motif in motifs.itertuples(index=False):
            pattern = str(motif.pattern).strip().upper().replace("T", "U")
            if not pattern or set(pattern) - IUPAC_RNA_ALPHABET:
                raise ValueError(f"unsupported motif alphabet: {motif.motif_id}")

            is_pure_acgu = set(pattern).issubset(set("ACGU"))
            p_len = len(pattern)

            if is_pure_acgu:
                for i in range(len(seq) - p_len + 1):
                    if seq[i:i + p_len] != pattern:
                        continue
                    region = coords[i:i + p_len]
                    if len(region) != p_len or sorted(region) != list(range(min(region), max(region) + 1)):
                        continue
                    rows.append({"regulatory_site_id": f"{iso.isoform_id}:{motif.motif_id}:{i}",
                                 "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                                 "motif_id": motif.motif_id, "regulator": motif.regulator,
                                 "regulator_class": motif.regulator_class,
                                 "isoform_start": i, "isoform_end": i + p_len,
                                 "chrom": pas.chrom, "genomic_start": min(region),
                                 "genomic_end": max(region), "strand": pas.strand,
                                 "evidence_tier": motif.evidence_tier, "source": motif.source,
                                 "version": motif.version, "source_url": motif.url,
                                 "pubmed_id": getattr(motif, "pubmed_id", ""),
                                 "assay": getattr(motif, "assay", ""),
                                 "cell_context": getattr(motif, "cell_context", ""),
                                 "matrix_id": getattr(motif, "matrix_id", "")})
            else:
                regex_str = "".join(IUPAC_TO_REGEX[c] for c in pattern)
                pattern_re = re.compile(f"(?=({regex_str}))")
                for m in pattern_re.finditer(seq):
                    i = m.start()
                    match_str = m.group(1)
                    match_len = len(match_str)
                    region = coords[i:i + match_len]
                    if len(region) != match_len or sorted(region) != list(range(min(region), max(region) + 1)):
                        continue
                    rows.append({"regulatory_site_id": f"{iso.isoform_id}:{motif.motif_id}:{i}",
                                 "isoform_id": iso.isoform_id, "site_id": iso.site_id,
                                 "motif_id": motif.motif_id, "regulator": motif.regulator,
                                 "regulator_class": motif.regulator_class,
                                 "isoform_start": i, "isoform_end": i + match_len,
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
