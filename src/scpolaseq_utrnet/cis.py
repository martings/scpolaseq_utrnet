"""Cis-regulatory motif dosage and dynamic survival field modeling across alternative 3'UTRs."""
from __future__ import annotations

import logging
from typing import Any
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def annotate_motif_survival_probabilities(
    regulatory_sites: pd.DataFrame,
    pas_sites: pd.DataFrame,
    profiles: pd.DataFrame,
    population_col: str = "celltype_corrected",
    usage_col: str = "relative_isoform_usage",
) -> pd.DataFrame:
    """Annotate each regulatory motif match with its continuous survival probability S(m; c).

    For strand '+': S(m; c) = sum_{k: x_k >= m_end} theta_k(c)
    For strand '-': S(m; c) = sum_{k: x_k <= m_start} theta_k(c)

    Parameters
    ----------
    regulatory_sites : pd.DataFrame
        Table of scanned regulatory sites with columns:
        ['gene_id', 'regulator', 'regulator_class', 'motif_id', 'chrom',
         'genomic_start', 'genomic_end', 'strand'].
    pas_sites : pd.DataFrame
        PAS reference sites with columns ['site_id', 'start', 'end', 'strand'].
    profiles : pd.DataFrame
        Population-specific PAS usage with columns ['site_id', population_col, usage_col].
    population_col : str
        Column denoting cell population/condition.
    usage_col : str
        Column denoting relative isoform usage theta_p(c) (sum to 1 per population).

    Returns
    -------
    pd.DataFrame
        Table of motifs x populations with columns:
        ['gene_id', 'regulator', 'regulator_class', 'motif_id', 'chrom',
         'genomic_start', 'genomic_end', 'strand', 'region_type',
         population_col, 'survival_probability']
    """
    if regulatory_sites.empty or pas_sites.empty or profiles.empty:
        return pd.DataFrame(columns=[
            "gene_id", "regulator", "regulator_class", "motif_id", "chrom",
            "genomic_start", "genomic_end", "strand", "region_type",
            population_col, "survival_probability"
        ])

    # Ensure pas_sites coordinates are numeric
    pas_df = pas_sites[["site_id", "start", "end", "strand"]].copy().drop_duplicates("site_id")
    pas_df["cleavage_coord"] = pas_df["start"].astype(int)

    # Merge PAS positions with population usage
    usage_df = profiles[[population_col, "site_id", usage_col]].copy()
    usage_df = usage_df.merge(pas_df[["site_id", "cleavage_coord", "strand"]], on="site_id", how="inner")

    # Determine proximal PAS per gene/strand
    strand = str(regulatory_sites["strand"].iloc[0])
    if strand == "-":
        proximal_coord = pas_df["cleavage_coord"].max()
        distal_coord = pas_df["cleavage_coord"].min()
    else:
        proximal_coord = pas_df["cleavage_coord"].min()
        distal_coord = pas_df["cleavage_coord"].max()

    # Pre-calculate survival field for each unique motif position
    unique_motifs = regulatory_sites[[
        "regulatory_site_id", "gene_id", "regulator", "regulator_class",
        "motif_id", "chrom", "genomic_start", "genomic_end", "strand"
    ]].drop_duplicates("regulatory_site_id").copy()

    unique_motifs["genomic_start"] = unique_motifs["genomic_start"].astype(int)
    unique_motifs["genomic_end"] = unique_motifs["genomic_end"].astype(int)

    # Classify region: cUTR (constitutive common region) vs aUTR (alternative dynamic region)
    if strand == "-":
        unique_motifs["region_type"] = np.where(
            unique_motifs["genomic_start"] >= proximal_coord, "cUTR",
            np.where(unique_motifs["genomic_end"] <= distal_coord, "post_distal", "aUTR")
        )
    else:
        unique_motifs["region_type"] = np.where(
            unique_motifs["genomic_end"] <= proximal_coord, "cUTR",
            np.where(unique_motifs["genomic_start"] >= distal_coord, "post_distal", "aUTR")
        )

    # Filter out post_distal motifs (they are beyond all known PAS sites)
    valid_motifs = unique_motifs[unique_motifs["region_type"] != "post_distal"].copy()
    if valid_motifs.empty:
        return pd.DataFrame()

    # Group usage by population for fast vectorization
    records = []
    for pop, pop_usage in usage_df.groupby(population_col):
        coords = pop_usage["cleavage_coord"].values
        thetas = pop_usage[usage_col].fillna(0.0).values
        total_theta = thetas.sum()
        norm_thetas = thetas / total_theta if total_theta > 0 else thetas

        for m in valid_motifs.itertuples(index=False):
            if m.region_type == "cUTR":
                surv = 1.0
            else:
                if strand == "-":
                    # For '-' strand, transcription flows high -> low. Retained if cleavage <= m.genomic_start
                    mask = coords <= m.genomic_start
                else:
                    # For '+' strand, transcription flows low -> high. Retained if cleavage >= m.genomic_end
                    mask = coords >= m.genomic_end
                surv = float(norm_thetas[mask].sum())

            records.append({
                "regulatory_site_id": m.regulatory_site_id,
                "gene_id": m.gene_id,
                "regulator": m.regulator,
                "regulator_class": m.regulator_class,
                "motif_id": m.motif_id,
                "chrom": m.chrom,
                "genomic_start": m.genomic_start,
                "genomic_end": m.genomic_end,
                "strand": m.strand,
                "region_type": m.region_type,
                population_col: pop,
                "survival_probability": min(max(surv, 0.0), 1.0),
            })

    return pd.DataFrame(records)


def compute_fractional_cis_dosage(
    motif_survival_df: pd.DataFrame,
    population_col: str = "celltype_corrected",
) -> pd.DataFrame:
    """Compute the fractional expected dosage <N_r(g, c)> of regulatory binding sites per mRNA molecule.

    Decomposes dosage into constitutive common region <N_r^{cUTR}> and dynamic alternative region <N_r^{aUTR}>.

    Parameters
    ----------
    motif_survival_df : pd.DataFrame
        Output of `annotate_motif_survival_probabilities`.
    population_col : str
        Column denoting cell population/condition.

    Returns
    -------
    pd.DataFrame
        Table aggregated per (gene_id, regulator, regulator_class, population):
        ['gene_id', 'regulator', 'regulator_class', population_col,
         'n_motifs_cutr', 'dosage_autr', 'total_dosage', 'autr_dosage_fraction']
    """
    if motif_survival_df.empty:
        return pd.DataFrame(columns=[
            "gene_id", "regulator", "regulator_class", population_col,
            "n_motifs_cutr", "dosage_autr", "total_dosage", "autr_dosage_fraction"
        ])

    rows = []
    group_cols = ["gene_id", "regulator", "regulator_class", population_col]

    for (gene_id, reg, reg_cls, pop), grp in motif_survival_df.groupby(group_cols):
        cutr_mask = grp["region_type"] == "cUTR"
        autr_mask = grp["region_type"] == "aUTR"

        # cUTR motifs are present with survival = 1.0
        n_cutr = int(cutr_mask.sum())

        # aUTR motifs have dynamic fractional survival S(m; c) in [0, 1]
        dosage_autr = float(grp.loc[autr_mask, "survival_probability"].sum())
        total_dosage = float(n_cutr + dosage_autr)
        frac = dosage_autr / total_dosage if total_dosage > 0 else 0.0

        rows.append({
            "gene_id": gene_id,
            "regulator": reg,
            "regulator_class": reg_cls,
            population_col: pop,
            "n_motifs_cutr": n_cutr,
            "dosage_autr": dosage_autr,
            "total_dosage": total_dosage,
            "autr_dosage_fraction": frac,
        })

    return pd.DataFrame(rows)


def compute_dynamic_dosage_flux(
    dosage_df: pd.DataFrame,
    start_population: str,
    end_population: str,
    population_col: str = "celltype_corrected",
    delta_t: float = 1.0,
) -> pd.DataFrame:
    """Compute step-wise dynamic flux / rate of change of cis-dosage: Delta <N_r> / Delta t.

    Parameters
    ----------
    dosage_df : pd.DataFrame
        Output of `compute_fractional_cis_dosage`.
    start_population : str
        Reference or baseline state (e.g., 'DF_day1').
    end_population : str
        Target or follow-up state (e.g., 'DF_wk2').
    population_col : str
        Column denoting population name.
    delta_t : float
        Elapsed time interval.

    Returns
    -------
    pd.DataFrame
        Table of dosage deltas:
        ['gene_id', 'regulator', 'regulator_class', 'start_state', 'end_state',
         'dosage_start', 'dosage_end', 'delta_dosage', 'dosage_velocity',
         'dosage_flux_class']
    """
    if dosage_df.empty:
        return pd.DataFrame(columns=[
            "gene_id", "regulator", "regulator_class", "start_state", "end_state",
            "dosage_start", "dosage_end", "delta_dosage", "dosage_velocity",
            "dosage_flux_class"
        ])

    sub_start = dosage_df[dosage_df[population_col] == start_population].set_index(
        ["gene_id", "regulator", "regulator_class"]
    )
    sub_end = dosage_df[dosage_df[population_col] == end_population].set_index(
        ["gene_id", "regulator", "regulator_class"]
    )

    common_keys = sub_start.index.intersection(sub_end.index)
    records = []

    for key in common_keys:
        gene_id, reg, reg_cls = key
        d0 = float(sub_start.loc[key, "total_dosage"])
        d1 = float(sub_end.loc[key, "total_dosage"])
        delta_d = d1 - d0
        vel = delta_d / delta_t if delta_t > 0 else delta_d

        if abs(delta_d) < 1e-4:
            flux_cls = "invariant_dosage"
        elif delta_d > 0.05:
            flux_cls = "cis_motif_expansion"
        elif delta_d < -0.05:
            flux_cls = "cis_motif_shedding"
        else:
            flux_cls = "minor_modulation"

        records.append({
            "gene_id": gene_id,
            "regulator": reg,
            "regulator_class": reg_cls,
            "start_state": start_population,
            "end_state": end_population,
            "dosage_start": d0,
            "dosage_end": d1,
            "delta_dosage": delta_d,
            "dosage_velocity": vel,
            "dosage_flux_class": flux_cls,
        })

    return pd.DataFrame(records)
