"""Strand-aware, ambiguity-preserving 3′UTR candidate reconstruction."""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import re
import pysam
import pandas as pd


ATTR = re.compile(r'(\w+) "([^"]*)"')
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_gene_models(gtf: Path, gene_id: str) -> dict[str, dict]:
    models: dict[str, dict] = defaultdict(lambda: {"exon": [], "CDS": [], "UTR": [], "meta": {}})
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#") or f'gene_id "{gene_id}' not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attr = dict(ATTR.findall(fields[8]))
            if attr.get("gene_id", "").split(".")[0] != gene_id or "transcript_id" not in attr:
                continue
            model = models[attr["transcript_id"]]
            model["meta"] = {"chrom": fields[0], "strand": fields[6],
                             "transcript_type": attr.get("transcript_type", ""),
                             "gene_id": attr["gene_id"]}
            if fields[2] in {"exon", "CDS", "UTR"}:
                model[fields[2]].append((int(fields[3]), int(fields[4])))
    return dict(models)


def _utr_intervals(model: dict) -> list[tuple[int, int]]:
    cds = model["CDS"]
    if not cds:
        return []
    strand = model["meta"]["strand"]
    boundary = max(end for _, end in cds) if strand == "+" else min(start for start, _ in cds)
    intervals = []
    for start, end in model["UTR"]:
        if strand == "+" and end > boundary:
            intervals.append((max(start, boundary + 1), end))
        elif strand == "-" and start < boundary:
            intervals.append((start, min(end, boundary - 1)))
    return sorted(set((s, e) for s, e in intervals if s <= e), reverse=strand == "-")


def _sequence_to_pas(fasta: pysam.FastaFile, chrom: str, strand: str,
                     intervals: list[tuple[int, int]], endpoint: int) -> tuple[str, list[int]]:
    pieces: list[str] = []
    coords: list[int] = []
    for start, end in intervals:
        if strand == "+":
            if start > endpoint:
                break
            end = min(end, endpoint)
            if start > end:
                continue
            pieces.append(fasta.fetch(chrom, start - 1, end).upper())
            coords.extend(range(start, end + 1))
        else:
            if end < endpoint:
                break
            start = max(start, endpoint)
            if start > end:
                continue
            pieces.append(fasta.fetch(chrom, start - 1, end).translate(COMPLEMENT)[::-1].upper())
            coords.extend(range(end, start - 1, -1))
    return "".join(pieces), coords


def reconstruct_isoforms(sites: pd.DataFrame, gtf: Path, fasta_path: Path,
                         max_span: int = 50000) -> tuple[pd.DataFrame, dict[str, list[int]]]:
    gene_id = str(sites.gene_id.iloc[0]).split(".")[0]
    models = parse_gene_models(gtf, gene_id)
    records = []
    coordinate_maps: dict[str, list[int]] = {}
    gene_pas_span = int(sites.start.max() - sites.start.min())
    span_status = "suspicious_gene_pas_span" if gene_pas_span > max_span else "within_span_limit"
    with pysam.FastaFile(str(fasta_path)) as fasta:
        for site in sites.itertuples(index=False):
            endpoint = int(site.start)
            candidates = []
            for transcript, model in models.items():
                meta = model["meta"]
                if meta["strand"] != site.strand or meta["chrom"] != site.chrom or not model["exon"]:
                    continue
                terminal = (max(model["exon"], key=lambda x: x[1]) if site.strand == "+"
                            else min(model["exon"], key=lambda x: x[0]))
                if not terminal[0] <= endpoint <= terminal[1]:
                    continue
                intervals = _utr_intervals(model)
                if not any(a <= endpoint <= b for a, b in intervals):
                    continue
                seq, coords = _sequence_to_pas(fasta, site.chrom, site.strand, intervals, endpoint)
                if seq and len(seq) == len(coords) and len(seq) <= max_span:
                    candidates.append((transcript, seq, coords, meta["transcript_type"], terminal))
            distinct = {}
            for transcript, seq, coords, kind, terminal in candidates:
                distinct.setdefault(seq, []).append((transcript, coords, kind, terminal))
            if not distinct:
                records.append({"isoform_id": f"utr:{site.pas_reference_id}:unresolved", "site_id": site.site_id,
                                "pas_reference_id": site.pas_reference_id, "gene_id": gene_id,
                                "transcript_ids": "", "transcript_types": "", "terminal_exon": "",
                                "sequence": "", "sequence_length": 0, "ambiguity_status": "no_compatible_3utr_model",
                                "model_weight": 0.0, "gene_pas_span_bp": gene_pas_span,
                                "span_status": span_status})
                continue
            for i, (seq, group) in enumerate(sorted(distinct.items(), key=lambda x: (-len(x[0]), x[0])), 1):
                iso_id = f"utr:{site.pas_reference_id}:{i}"
                status = "unique" if len(candidates) == 1 else (
                    "same_sequence_multiple_transcripts" if len(distinct) == 1 else "ambiguous_sequence_models")
                records.append({"isoform_id": iso_id, "site_id": site.site_id,
                                "pas_reference_id": site.pas_reference_id, "gene_id": gene_id,
                                "transcript_ids": ";".join(x[0] for x in group),
                                "transcript_types": ";".join(sorted(set(x[2] for x in group))),
                                "terminal_exon": f"{site.chrom}:{group[0][3][0]}-{group[0][3][1]}:{site.strand}",
                                "sequence": seq, "sequence_length": len(seq),
                                "ambiguity_status": status, "model_weight": 1.0 / len(distinct),
                                "gene_pas_span_bp": gene_pas_span, "span_status": span_status})
                coordinate_maps[iso_id] = group[0][1]
    return pd.DataFrame.from_records(records), coordinate_maps
