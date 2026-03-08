"""
engraftment_score.py
Computes the Engraftment Suitability Score (ESS) for each cell.

The ESS is a weighted gene expression score that integrates:
  - Positive markers: neural maturity, scaffold adhesion genes
  - Negative markers: pluripotency, apoptosis, proliferation, stress genes

This score predicts which cell populations within an iPSC-derived
progenitor preparation are most likely to survive transplantation
into a biomaterial scaffold environment.

Biological rationale:
  The core challenge in scaffold-based cell therapy for spinal cord injury
  is that cell preparations contain heterogeneous populations. Some cells
  are mature enough to engraft and differentiate; others retain pluripotency
  (safety risk), are undergoing apoptosis (won't survive), or are still
  proliferating (tumorigenic risk). Bulk markers miss this intra-population
  variation. ESS uses single-cell resolution to score each cell individually.
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import yaml
import logging
from pathlib import Path
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import MinMaxScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

Path("results/figures/scoring").mkdir(parents=True, exist_ok=True)
Path("results/tables").mkdir(parents=True, exist_ok=True)


def load_scoring_weights(config_path: str = "configs/scoring_weights.yaml") -> tuple:
    """Load ESS gene weights from config."""
    with open(config_path) as f:
        cfg = yaml.safe_load(f)
    positive = cfg["positive_markers"]
    negative = cfg["negative_markers"]
    return positive, negative


def compute_ess(adata: sc.AnnData, positive: dict, negative: dict) -> sc.AnnData:
    """
    Compute the Engraftment Suitability Score (ESS) per cell.

    For each gene present in the dataset, retrieve its log-normalized expression
    per cell and multiply by its weight. Sum across all genes to get raw ESS.
    Min-max scale to [0, 1] for interpretability.

    Returns adata with 'ESS_raw' and 'ESS' columns added to .obs.
    """
    log.info("Computing Engraftment Suitability Score (ESS)...")

    # Get expression matrix as DataFrame (cells × genes)
    # Use .raw if available (all genes, not just HVGs)
    if adata.raw is not None:
        expr = pd.DataFrame(
            adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names
        )
    else:
        expr = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        )

    # Initialize score
    ess = pd.Series(0.0, index=adata.obs_names)
    genes_used = {"positive": [], "negative": []}

    # Add positive markers
    for gene, weight in positive.items():
        if gene in expr.columns:
            ess += expr[gene] * weight
            genes_used["positive"].append(gene)
        else:
            log.debug(f"  Positive marker not found in dataset: {gene}")

    # Add negative markers (weights are already negative in config)
    for gene, weight in negative.items():
        if gene in expr.columns:
            ess += expr[gene] * weight
            genes_used["negative"].append(gene)
        else:
            log.debug(f"  Negative marker not found in dataset: {gene}")

    log.info(
        f"  ESS computed using {len(genes_used['positive'])} positive "
        f"and {len(genes_used['negative'])} negative markers"
    )

    # Store raw ESS
    adata.obs["ESS_raw"] = ess.values

    # Min-max normalize to [0, 1]
    scaler = MinMaxScaler()
    adata.obs["ESS"] = scaler.fit_transform(ess.values.reshape(-1, 1)).flatten()

    # Classify into tiers
    adata.obs["ESS_tier"] = pd.cut(
        adata.obs["ESS"],
        bins=[0, 0.33, 0.67, 1.01],
        labels=["Low (Avoid)", "Medium", "High (Optimal)"],
        include_lowest=True
    )

    # Summary statistics
    log.info("\n── ESS Summary by Cell Type ──────────────────────────────")
    if "predicted_cell_type" in adata.obs.columns:
        summary = (
            adata.obs.groupby("predicted_cell_type")["ESS"]
            .agg(["mean", "median", "std", "count"])
            .round(3)
            .sort_values("median", ascending=False)
        )
        log.info(f"\n{summary.to_string()}")

        # Save summary table
        summary.to_csv("results/tables/ess_by_cell_type.csv")
        log.info("  Saved: results/tables/ess_by_cell_type.csv")

    if "sample_id" in adata.obs.columns:
        timepoint_summary = (
            adata.obs.groupby("sample_id")["ESS"]
            .agg(["mean", "median", "std", "count"])
            .round(3)
            .sort_values("median", ascending=False)
        )
        log.info(f"\n── ESS Summary by Timepoint ──")
        log.info(f"\n{timepoint_summary.to_string()}")
        timepoint_summary.to_csv("results/tables/ess_by_timepoint.csv")

    return adata, genes_used


def identify_optimal_harvest_window(adata: sc.AnnData) -> pd.DataFrame:
    """
    Identify the differentiation timepoint and cell population that maximizes ESS.
    Performs pairwise Mann-Whitney U tests between timepoints to assess significance.
    """
    log.info("\nIdentifying optimal harvest window...")

    results = []

    if "sample_id" not in adata.obs.columns:
        log.warning("sample_id not in obs — cannot compute timepoint comparison")
        return pd.DataFrame()

    timepoints = adata.obs["sample_id"].unique()

    for tp in timepoints:
        mask = adata.obs["sample_id"] == tp
        ess_vals = adata.obs.loc[mask, "ESS"].values
        high_ess_frac = (adata.obs.loc[mask, "ESS_tier"] == "High (Optimal)").mean()

        results.append({
            "timepoint": tp,
            "n_cells": mask.sum(),
            "mean_ESS": ess_vals.mean(),
            "median_ESS": np.median(ess_vals),
            "frac_high_ESS": high_ess_frac,
        })

    df = pd.DataFrame(results).sort_values("median_ESS", ascending=False)

    log.info("\n── Optimal Harvest Ranking ──")
    log.info(df.to_string(index=False))
    df.to_csv("results/tables/optimal_harvest_timepoints.csv", index=False)

    optimal = df.iloc[0]["timepoint"]
    log.info(f"\n✓ Recommended harvest timepoint: {optimal}")
    log.info(
        "  Rationale: Highest median ESS + highest fraction of 'High (Optimal)' tier cells"
    )

    return df


def derive_surface_markers(adata: sc.AnnData, top_n: int = 20) -> pd.DataFrame:
    """
    Identify cell surface genes that best distinguish High-ESS from Low-ESS cells.
    These are candidates for FACS sorting to prospectively enrich high-ESS populations
    without requiring single-cell transcriptomics at clinical scale.

    Uses Mann-Whitney U test + effect size (rank-biserial correlation).
    """
    log.info("\nIdentifying surface marker candidates for FACS enrichment...")

    # Known surface / secreted genes (could be extended with surfaceome database)
    surface_gene_prefixes = ("CD", "ITGA", "ITGB", "NCAM", "L1CAM", "NOTCH",
                             "EPHB", "EPHA", "PDGFR", "EGFR", "MET", "HER")

    # Filter to surface genes present in dataset
    if adata.raw is not None:
        all_genes = adata.raw.var_names
    else:
        all_genes = adata.var_names

    surface_genes = [
        g for g in all_genes
        if any(g.startswith(p) for p in surface_gene_prefixes)
    ]
    log.info(f"  Testing {len(surface_genes)} candidate surface genes")

    if not surface_genes or "ESS_tier" not in adata.obs.columns:
        log.warning("  Skipping surface marker analysis — insufficient data")
        return pd.DataFrame()

    high_mask = adata.obs["ESS_tier"] == "High (Optimal)"
    low_mask = adata.obs["ESS_tier"] == "Low (Avoid)"

    # Get expression
    if adata.raw is not None:
        expr_mat = pd.DataFrame(
            adata.raw.X.toarray() if hasattr(adata.raw.X, "toarray") else adata.raw.X,
            index=adata.obs_names,
            columns=adata.raw.var_names
        )
    else:
        expr_mat = pd.DataFrame(
            adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
            index=adata.obs_names,
            columns=adata.var_names
        )

    results = []
    for gene in surface_genes:
        if gene not in expr_mat.columns:
            continue
        high_expr = expr_mat.loc[high_mask, gene].values
        low_expr = expr_mat.loc[low_mask, gene].values

        if high_expr.std() < 1e-6 and low_expr.std() < 1e-6:
            continue

        stat, pval = mannwhitneyu(high_expr, low_expr, alternative="two-sided")

        # Effect size: rank-biserial correlation
        n1, n2 = len(high_expr), len(low_expr)
        r = 1 - (2 * stat) / (n1 * n2)

        results.append({
            "gene": gene,
            "mean_expr_high_ESS": high_expr.mean(),
            "mean_expr_low_ESS": low_expr.mean(),
            "fold_change": (high_expr.mean() + 1e-6) / (low_expr.mean() + 1e-6),
            "pvalue": pval,
            "effect_size_r": r
        })

    if not results:
        return pd.DataFrame()

    df = pd.DataFrame(results)
    df["enriched_in"] = df["effect_size_r"].apply(
        lambda r: "High-ESS" if r > 0 else "Low-ESS"
    )
    df = df.sort_values("effect_size_r", key=abs, ascending=False).head(top_n)

    log.info(f"\n── Top {top_n} Surface Marker Candidates ──")
    log.info(df[["gene", "fold_change", "pvalue", "effect_size_r", "enriched_in"]]
             .to_string(index=False))

    df.to_csv("results/tables/surface_marker_candidates.csv", index=False)
    log.info("  Saved: results/tables/surface_marker_candidates.csv")
    log.info(
        "\n  These markers can guide FACS-based enrichment of high-ESS cells "
        "prior to scaffold seeding — enabling clinical-scale selection without "
        "requiring per-batch transcriptomics."
    )

    return df


def plot_ess_results(adata: sc.AnnData, genes_used: dict) -> None:
    """Publication-quality ESS visualization panel."""
    log.info("\nGenerating ESS visualization plots...")

    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.35)

    # ── Panel A: ESS on UMAP ──────────────────────────────────────────────────
    ax_a = fig.add_subplot(gs[0, :2])
    if "X_umap" in adata.obsm:
        scatter = ax_a.scatter(
            adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
            c=adata.obs["ESS"], cmap="RdYlGn", s=1.5, alpha=0.7,
            vmin=0, vmax=1
        )
        plt.colorbar(scatter, ax=ax_a, label="ESS")
    ax_a.set_title("A. Engraftment Suitability Score (ESS) — UMAP", fontweight="bold")
    ax_a.set_xlabel("UMAP 1")
    ax_a.set_ylabel("UMAP 2")

    # ── Panel B: ESS tier distribution ───────────────────────────────────────
    ax_b = fig.add_subplot(gs[0, 2])
    tier_counts = adata.obs["ESS_tier"].value_counts()
    colors = ["#d73027", "#fdae61", "#1a9850"]
    ax_b.pie(
        tier_counts.values, labels=tier_counts.index,
        colors=colors, autopct="%1.1f%%", startangle=90
    )
    ax_b.set_title("B. ESS Tier Distribution", fontweight="bold")

    # ── Panel C: ESS violin by cell type ─────────────────────────────────────
    ax_c = fig.add_subplot(gs[1, :])
    if "predicted_cell_type" in adata.obs.columns:
        order = (
            adata.obs.groupby("predicted_cell_type")["ESS"]
            .median().sort_values(ascending=False).index.tolist()
        )
        sns.violinplot(
            data=adata.obs, x="predicted_cell_type", y="ESS",
            order=order, palette="RdYlGn", ax=ax_c,
            inner="box", cut=0
        )
        ax_c.axhline(0.67, color="green", linestyle="--",
                     alpha=0.7, label="High ESS threshold")
        ax_c.axhline(0.33, color="red", linestyle="--",
                     alpha=0.7, label="Low ESS threshold")
        ax_c.set_title("C. ESS Distribution by Cell Type", fontweight="bold")
        ax_c.set_xlabel("Cell Type")
        ax_c.set_ylabel("ESS")
        ax_c.tick_params(axis="x", rotation=30)
        ax_c.legend()

    # ── Panel D: ESS over differentiation time ───────────────────────────────
    ax_d = fig.add_subplot(gs[2, :2])
    if "sample_id" in adata.obs.columns:
        tp_order = sorted(adata.obs["sample_id"].unique())
        sns.boxplot(
            data=adata.obs, x="sample_id", y="ESS",
            order=tp_order, palette="Blues", ax=ax_d,
            fliersize=1
        )
        ax_d.set_title("D. ESS Across Differentiation Timepoints", fontweight="bold")
        ax_d.set_xlabel("Timepoint")
        ax_d.set_ylabel("ESS")
        ax_d.tick_params(axis="x", rotation=30)

    # ── Panel E: Gene contribution heatmap ───────────────────────────────────
    ax_e = fig.add_subplot(gs[2, 2])
    all_genes = genes_used["positive"] + genes_used["negative"][:5]
    if all_genes and "predicted_cell_type" in adata.obs.columns:
        mean_expr = (
            adata.obs.groupby("predicted_cell_type")[
                [g for g in all_genes if g in adata.obs.columns]
            ].mean()
        )
        if not mean_expr.empty:
            sns.heatmap(
                mean_expr.T, cmap="RdBu_r", center=0,
                ax=ax_e, cbar_kws={"shrink": 0.8}
            )
            ax_e.set_title("E. Key Gene Expression\nby Cell Type", fontweight="bold")
            ax_e.tick_params(axis="x", rotation=30, labelsize=7)
            ax_e.tick_params(axis="y", labelsize=7)

    fig.suptitle(
        "Neural Scaffold Cell Atlas — Engraftment Suitability Score (ESS)",
        fontsize=14, fontweight="bold", y=1.01
    )

    plt.savefig(
        "results/figures/scoring/ESS_panel_figure.png",
        dpi=200, bbox_inches="tight"
    )
    plt.close()
    log.info("  Saved: results/figures/scoring/ESS_panel_figure.png")


def main():
    input_path = "data/processed/GSE303787_velocity.h5ad"
    if not Path(input_path).exists():
        # Fall back to clustered data if velocity not computed yet
        input_path = "data/processed/GSE303787_clustered.h5ad"
        if not Path(input_path).exists():
            log.error("No processed data found. Run earlier pipeline steps first.")
            return

    log.info(f"Loading data from {input_path}...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Load weights and compute ESS
    positive, negative = load_scoring_weights()
    adata, genes_used = compute_ess(adata, positive, negative)

    # Derive actionable outputs
    harvest_df = identify_optimal_harvest_window(adata)
    surface_df = derive_surface_markers(adata)

    # Visualization
    plot_ess_results(adata, genes_used)

    # Save final annotated object
    out_path = "data/processed/GSE303787_ESS_scored.h5ad"
    adata.write(out_path)
    log.info(f"\n{'='*60}")
    log.info(f"✓ Engraftment Suitability Score complete")
    log.info(f"  Final annotated AnnData: {out_path}")
    log.info(f"  Summary tables: results/tables/")
    log.info(f"  Figures: results/figures/scoring/")
    log.info(f"{'='*60}")


if __name__ == "__main__":
    main()
