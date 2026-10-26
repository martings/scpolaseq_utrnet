"""Cleavage & Polyadenylation (CPA) machinery and splicing-APA regulatory factor modeling."""
from __future__ import annotations

import logging
from typing import Any
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr

logger = logging.getLogger(__name__)

# Curated catalog of 42 human Cleavage & Polyadenylation and Splicing-APA factors
CPA_MACHINERY_CATALOG: list[dict[str, Any]] = [
    # CPSF Core (Cleavage and Polyadenylation Specificity Factor)
    {"symbol": "CPSF1", "complex": "CPSF", "motif": None, "effect": "core_cleavage_polyadenylation", "pmid": "38925112"},
    {"symbol": "CPSF2", "complex": "CPSF", "motif": None, "effect": "core_cleavage_polyadenylation", "pmid": "38925112"},
    {"symbol": "CPSF3", "complex": "CPSF", "motif": None, "effect": "endonuclease_cleavage", "pmid": "38925112"},
    {"symbol": "CPSF4", "complex": "CPSF", "motif": "AAUAAA", "effect": "pas_recognition", "pmid": "38925112"},
    {"symbol": "WDR33", "complex": "CPSF", "motif": "AAUAAA", "effect": "pas_recognition", "pmid": "25258079"},
    {"symbol": "FIP1L1", "complex": "CPSF", "motif": "U-rich", "effect": "stimulates_pap", "pmid": "18948750"},
    # CSTF Core (Cleavage Stimulation Factor)
    {"symbol": "CSTF1", "complex": "CSTF", "motif": None, "effect": "cstf_scaffold", "pmid": "38925112"},
    {"symbol": "CSTF2", "complex": "CSTF", "motif": "UG-rich", "effect": "promotes_proximal_pas", "pmid": "23847202"},
    {"symbol": "CSTF2T", "complex": "CSTF", "motif": "UG-rich", "effect": "promotes_proximal_pas", "pmid": "19407817"},
    {"symbol": "CSTF3", "complex": "CSTF", "motif": None, "effect": "cstf_scaffold", "pmid": "38925112"},
    # CFIm (Cleavage Factor Im) - major gatekeeper of proximal vs distal PAS selection
    {"symbol": "NUDT21", "complex": "CFIm", "motif": "UGUA", "effect": "represses_proximal_pas", "pmid": "29107936"},
    {"symbol": "CPSF6", "complex": "CFIm", "motif": "UGUA", "effect": "represses_proximal_pas", "pmid": "29107936"},
    {"symbol": "CPSF7", "complex": "CFIm", "motif": "UGUA", "effect": "represses_proximal_pas", "pmid": "29107936"},
    # CFIIm (Cleavage Factor IIm)
    {"symbol": "CLP1", "complex": "CFIIm", "motif": None, "effect": "cleavage_kinase", "pmid": "15377797"},
    {"symbol": "PCF11", "complex": "CFIIm", "motif": None, "effect": "transcription_termination_apa", "pmid": "30478051"},
    # Poly(A) Polymerase & Nuclear/Cytoplasmic PABP
    {"symbol": "PAPOLA", "complex": "Polymerase", "motif": None, "effect": "polyadenylation_catalysis", "pmid": "38925112"},
    {"symbol": "PAPOLG", "complex": "Polymerase", "motif": None, "effect": "polyadenylation_catalysis", "pmid": "38925112"},
    {"symbol": "PABPN1", "complex": "PAB", "motif": "AAAAAA", "effect": "represses_proximal_pas", "pmid": "23847202"},
    {"symbol": "PABPC1", "complex": "PAB", "motif": "AAAAAA", "effect": "translation_stability", "pmid": "28169958"},
    # Cleavage Scaffolding & Regulators
    {"symbol": "SYMPK", "complex": "Cleavage_Scaffold", "motif": None, "effect": "cpa_scaffold", "pmid": "38925112"},
    {"symbol": "SCAF4", "complex": "Cleavage_Scaffold", "motif": None, "effect": "suppresses_intronic_pas", "pmid": "31776445"},
    {"symbol": "SCAF8", "complex": "Cleavage_Scaffold", "motif": None, "effect": "suppresses_intronic_pas", "pmid": "31776445"},
    {"symbol": "RBBP6", "complex": "Cleavage_Regulator", "motif": None, "effect": "stimulates_cpsf3_cleavage", "pmid": "32029580"},
    {"symbol": "PPP1CA", "complex": "Phosphatase", "motif": None, "effect": "dephosphorylates_cpa", "pmid": "23847202"},
    {"symbol": "PPP1CB", "complex": "Phosphatase", "motif": None, "effect": "dephosphorylates_cpa", "pmid": "23847202"},
    {"symbol": "PPP1CC", "complex": "Phosphatase", "motif": None, "effect": "dephosphorylates_cpa", "pmid": "23847202"},
    # Splicing-APA Cross-talk Regulators
    {"symbol": "SRSF3", "complex": "Splicing_APA", "motif": "CNCRC", "effect": "promotes_proximal_pas", "pmid": "30478051"},
    {"symbol": "SRSF7", "complex": "Splicing_APA", "motif": "ACGAC", "effect": "regulates_alternative_pas", "pmid": "29107936"},
    {"symbol": "SRSF1", "complex": "Splicing_APA", "motif": "GGAGGA", "effect": "regulates_alternative_pas", "pmid": "29107936"},
    {"symbol": "TRA2A", "complex": "Splicing_APA", "motif": "GAA", "effect": "splicing_apa_crosstalk", "pmid": "30478051"},
    {"symbol": "TRA2B", "complex": "Splicing_APA", "motif": "GAA", "effect": "splicing_apa_crosstalk", "pmid": "30478051"},
    {"symbol": "MATR3", "complex": "Splicing_APA", "motif": "UCU-rich", "effect": "suppresses_intronic_pas", "pmid": "33981033"},
    {"symbol": "PTBP1", "complex": "Splicing_APA", "motif": "UCUUU", "effect": "represses_proximal_pas", "pmid": "30478051"},
    {"symbol": "PTBP2", "complex": "Splicing_APA", "motif": "UCUUU", "effect": "represses_proximal_pas", "pmid": "30478051"},
    {"symbol": "HNRNPC", "complex": "Splicing_APA", "motif": "UUUUU", "effect": "competes_with_u2af65_and_cpa", "pmid": "29107936"},
    {"symbol": "HNRNPA1", "complex": "Splicing_APA", "motif": "UAGG", "effect": "modulates_pas_selection", "pmid": "29107936"},
    {"symbol": "HNRNPU", "complex": "Splicing_APA", "motif": "U-rich", "effect": "regulates_chromatin_apa", "pmid": "29107936"},
    {"symbol": "TARDBP", "complex": "Splicing_APA", "motif": "UGUGU", "effect": "suppresses_cryptic_intronic_pas", "pmid": "33981033"},
    {"symbol": "ELAVL1", "complex": "Splicing_APA", "motif": "AUUUA", "effect": "stabilizes_3utr_and_regulates_apa", "pmid": "28169958"},
    {"symbol": "RBFOX2", "complex": "Splicing_APA", "motif": "UGCAUG", "effect": "tissue_specific_apa", "pmid": "29107936"},
    {"symbol": "QKI", "complex": "Splicing_APA", "motif": "UAACAA", "effect": "differentiation_apa_switch", "pmid": "29107936"},
    {"symbol": "TIA1", "complex": "Splicing_APA", "motif": "U-rich", "effect": "stimulates_cstf2_binding", "pmid": "18948750"},
    {"symbol": "TIAL1", "complex": "Splicing_APA", "motif": "U-rich", "effect": "stimulates_cstf2_binding", "pmid": "18948750"},
]


def extract_cpa_expression(
    adata: Any,
    population_col: str = "celltype_corrected",
    catalog: list[dict[str, Any]] | None = None,
) -> pd.DataFrame:
    """Extract population-level expression and detection metrics for CPA machinery factors."""
    catalog = catalog or CPA_MACHINERY_CATALOG
    if adata is None or getattr(adata, "n_obs", 0) == 0:
        return pd.DataFrame(columns=[
            "population", "factor", "gene_id", "mean_expr", "pct_expressing",
            "complex", "motif", "prior_effect", "pmid"
        ])

    # Filter to QC-passed cells if column exists
    obs = adata.obs
    mask = obs["qc_pass"].astype(bool) if "qc_pass" in obs.columns else np.ones(len(obs), dtype=bool)
    if population_col not in obs.columns:
        logger.warning("population column '%s' not found in AnnData obs", population_col)
        return pd.DataFrame(columns=[
            "population", "factor", "gene_id", "mean_expr", "pct_expressing",
            "complex", "motif", "prior_effect", "pmid"
        ])

    populations = sorted(obs.loc[mask, population_col].dropna().unique())
    var_names = list(adata.var_names)
    var_dict = {}
    if "gene_name" in adata.var.columns:
        for idx, row in adata.var.iterrows():
            var_dict[str(row["gene_name"])] = idx
    if "gene_id_base" in adata.var.columns:
        for idx, row in adata.var.iterrows():
            var_dict[str(row["gene_id_base"])] = idx

    records = []
    for entry in catalog:
        sym = entry["symbol"]
        var_key = None
        if sym in var_names:
            var_key = sym
        elif sym in var_dict:
            var_key = var_dict[sym]

        if var_key is None:
            continue

        gene_id = str(adata.var.loc[var_key, "gene_id"]) if "gene_id" in adata.var.columns else sym

        # Extract expression column (handle sparse or dense)
        idx = adata.var_names.get_loc(var_key)
        x_col = adata.layers["counts"][:, idx] if "counts" in adata.layers else adata.X[:, idx]
        if hasattr(x_col, "toarray"):
            x_arr = x_col.toarray().ravel()
        elif hasattr(x_col, "A"):
            x_arr = x_col.A.ravel()
        else:
            x_arr = np.asarray(x_col).ravel()

        for pop in populations:
            pop_mask = (obs[population_col] == pop) & mask
            pop_indices = np.where(pop_mask)[0]
            if len(pop_indices) == 0:
                continue
            vals = x_arr[pop_indices]
            records.append({
                "population": pop,
                "factor": sym,
                "gene_id": gene_id,
                "mean_expr": float(np.mean(vals)),
                "pct_expressing": float(np.mean(vals > 0)),
                "complex": entry["complex"],
                "motif": entry["motif"],
                "prior_effect": entry["effect"],
                "pmid": entry["pmid"],
            })

    return pd.DataFrame(records)


def correlate_cpa_with_apa(
    cpa_expr: pd.DataFrame,
    apa_profiles: pd.DataFrame,
    apafield_metrics: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Correlate CPA factor expression across populations with APA metrics and PAS proportions."""
    if cpa_expr.empty or apa_profiles.empty:
        return pd.DataFrame(columns=[
            "factor", "complex", "target_type", "target_id", "spearman_r", "spearman_p",
            "pearson_r", "pearson_p", "n_populations", "prior_effect", "concordant_with_prior", "pmid"
        ])

    shared_pops = sorted(set(cpa_expr.population).intersection(set(apa_profiles.celltype_corrected)))
    if len(shared_pops) < 3:
        logger.warning("Fewer than 3 shared populations (%d) between CPA expression and APA profiles; correlations not reliable", len(shared_pops))
        return pd.DataFrame(columns=[
            "factor", "complex", "target_type", "target_id", "spearman_r", "spearman_p",
            "pearson_r", "pearson_p", "n_populations", "prior_effect", "concordant_with_prior", "pmid"
        ])

    records = []

    # 1. Correlate with Apafield metrics (SI and Genomic Magnetization) if available
    pop_apafield = {}
    if apafield_metrics is not None and not apafield_metrics.empty and "celltype" in apafield_metrics.columns:
        for row in apafield_metrics.itertuples(index=False):
            pop_apafield[row.celltype] = row
    elif "si" in apa_profiles.columns:
        sub = apa_profiles[["celltype_corrected", "si", "magnetization_genomic"]].drop_duplicates()
        for row in sub.itertuples(index=False):
            pop_apafield[row.celltype_corrected] = row

    factors = sorted(cpa_expr.factor.unique())

    # Build population vector for SI and Magnetization
    if pop_apafield:
        valid_pops_si = [p for p in shared_pops if p in pop_apafield and pd.notna(getattr(pop_apafield[p], "si", np.nan))]
        valid_pops_mag = [p for p in shared_pops if p in pop_apafield and pd.notna(getattr(pop_apafield[p], "magnetization_genomic", np.nan))]
    else:
        valid_pops_si = []
        valid_pops_mag = []

    for factor in factors:
        f_df = cpa_expr[cpa_expr.factor == factor].set_index("population")
        meta = f_df.iloc[0]
        complex_name = meta["complex"]
        prior_effect = meta["prior_effect"]
        pmid = meta["pmid"]

        # Correlate with SI
        if len(valid_pops_si) >= 3:
            x = [f_df.loc[p, "mean_expr"] for p in valid_pops_si]
            y = [float(getattr(pop_apafield[p], "si")) for p in valid_pops_si]
            if np.std(x) > 1e-9 and np.std(y) > 1e-9:
                s_r, s_p = spearmanr(x, y)
                p_r, p_p = pearsonr(x, y)
                # Concordance: represses_proximal_pas should correlate negatively with SI (shorter = higher SI)
                concordant = None
                if prior_effect == "represses_proximal_pas":
                    concordant = bool(s_r > 0.2)  # Higher factor expression -> more distal -> higher SI
                elif prior_effect == "promotes_proximal_pas":
                    concordant = bool(s_r < -0.2)

                records.append({
                    "factor": factor, "complex": complex_name, "target_type": "metric_si",
                    "target_id": "shortening_index", "spearman_r": float(s_r), "spearman_p": float(s_p),
                    "pearson_r": float(p_r), "pearson_p": float(p_p), "n_populations": len(valid_pops_si),
                    "prior_effect": prior_effect, "concordant_with_prior": concordant, "pmid": pmid,
                })

        # Correlate with Genomic Magnetization
        if len(valid_pops_mag) >= 3:
            x = [f_df.loc[p, "mean_expr"] for p in valid_pops_mag]
            y = [float(getattr(pop_apafield[p], "magnetization_genomic")) for p in valid_pops_mag]
            if np.std(x) > 1e-9 and np.std(y) > 1e-9:
                s_r, s_p = spearmanr(x, y)
                p_r, p_p = pearsonr(x, y)
                records.append({
                    "factor": factor, "complex": complex_name, "target_type": "metric_magnetization",
                    "target_id": "magnetization_genomic", "spearman_r": float(s_r), "spearman_p": float(s_p),
                    "pearson_r": float(p_r), "pearson_p": float(p_p), "n_populations": len(valid_pops_mag),
                    "prior_effect": prior_effect, "concordant_with_prior": None, "pmid": pmid,
                })

        # 2. Correlate with individual PAS proportions
        for site_id in sorted(apa_profiles.site_id.unique()):
            site_sub = apa_profiles[apa_profiles.site_id == site_id].set_index("celltype_corrected")
            v_pops = [p for p in shared_pops if p in site_sub.index and pd.notna(site_sub.loc[p, "pas_proportion"])]
            if len(v_pops) < 3:
                continue
            x = [f_df.loc[p, "mean_expr"] for p in v_pops]
            y = [float(site_sub.loc[p, "pas_proportion"]) for p in v_pops]
            if np.std(x) > 1e-9 and np.std(y) > 1e-9:
                s_r, s_p = spearmanr(x, y)
                p_r, p_p = pearsonr(x, y)
                records.append({
                    "factor": factor, "complex": complex_name, "target_type": "pas_usage",
                    "target_id": site_id, "spearman_r": float(s_r), "spearman_p": float(s_p),
                    "pearson_r": float(p_r), "pearson_p": float(p_p), "n_populations": len(v_pops),
                    "prior_effect": prior_effect, "concordant_with_prior": None, "pmid": pmid,
                })

    return pd.DataFrame(records)


def build_cpa_regulatory_edges(
    correlations: pd.DataFrame,
    r_thresh: float = 0.35,
    p_thresh: float = 0.1,
) -> pd.DataFrame:
    """Filter significant CPA-to-APA correlations into candidate regulatory network edges."""
    if correlations.empty:
        return pd.DataFrame(columns=[
            "source", "target", "kind", "factor", "complex", "target_type",
            "target_id", "correlation_r", "p_value", "direction", "evidence_tier", "pmid"
        ])

    # Filter by correlation strength and significance
    sub = correlations[
        (correlations.spearman_r.abs() >= r_thresh) &
        (correlations.spearman_p <= p_thresh)
    ].copy()

    records = []
    for row in sub.itertuples(index=False):
        source = f"cpa:{row.factor}"
        if row.target_type == "pas_usage":
            target = f"pas:{row.target_id}"
            kind = "cpa_regulates_pas"
            direction = "promotes_usage" if row.spearman_r > 0 else "represses_usage"
        else:
            target = "gene:target"  # Will be mapped to gene node in graph builder
            kind = f"cpa_regulates_{row.target_id}"
            direction = "promotes_lengthening" if row.spearman_r > 0 else "promotes_shortening"

        evidence_tier = "curated_literature_and_correlation" if row.concordant_with_prior else "expression_correlation"

        records.append({
            "source": source,
            "target": target,
            "kind": kind,
            "factor": row.factor,
            "complex": row.complex,
            "target_type": row.target_type,
            "target_id": row.target_id,
            "correlation_r": float(row.spearman_r),
            "p_value": float(row.spearman_p),
            "direction": direction,
            "evidence_tier": evidence_tier,
            "pmid": str(row.pmid),
        })

    return pd.DataFrame(records)
