import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

scv.settings.verbosity = 1
Path("results/figures/velocity").mkdir(parents=True, exist_ok=True)


def main():
    input_path = "data/processed/GSE303787_clustered.h5ad"
    if not Path(input_path).exists():
        log.error(f"Input not found: {input_path}. Run normalize_cluster.py first.")
        return

    log.info("Loading clustered data...")
    adata = sc.read_h5ad(input_path)
    log.info(f"  {adata.n_obs} cells x {adata.n_vars} genes")

    # ── Subsample to 10k cells to reduce RAM usage ────────────────────────────
    if adata.n_obs > 10000:
        log.info("Subsampling to 10,000 cells for velocity analysis...")
        sc.pp.subsample(adata, n_obs=10000, random_state=42)
        log.info(f"  Subsampled: {adata.n_obs} cells")

    # ── Simulate spliced/unspliced for pipeline demo ──────────────────────────
    # In real use: re-align with STARsolo --soloFeatures Gene Velocyto
    log.info("Setting up spliced/unspliced layers...")
    counts = adata.X.toarray() if hasattr(adata.X, "toarray") else np.array(adata.X)
    rng = np.random.default_rng(42)
    adata.layers["spliced"] = counts
    adata.layers["unspliced"] = np.abs(rng.normal(0, 0.1, counts.shape)) * counts * 0.3

    # ── scVelo pipeline ───────────────────────────────────────────────────────
    log.info("Filtering and normalizing for velocity...")
    scv.pp.filter_and_normalize(
        adata, min_shared_counts=20, enforce=True
    )

    log.info("Computing moments...")
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

    log.info("Fitting stochastic model (faster than dynamical)...")
    scv.tl.velocity(adata, mode="stochastic")

    log.info("Recomputing neighbors on subsampled data...")
    scv.pp.neighbors(adata, n_pcs=30, n_neighbors=30)

    log.info("Building velocity graph...")
    scv.tl.velocity_graph(adata, n_jobs=1)

    log.info("Computing UMAP for subsampled data...")
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_pca_harmony", random_state=42)
    sc.tl.umap(adata, random_state=42)

    # ── Plots ─────────────────────────────────────────────────────────────────
    log.info("Plotting velocity streams...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    scv.pl.velocity_embedding_stream(
        adata, basis="umap", color="sample_id",
        ax=axes[0], show=False, title="RNA Velocity — Timepoints"
    )
    if "predicted_cell_type" in adata.obs.columns:
        scv.pl.velocity_embedding_stream(
            adata, basis="umap", color="predicted_cell_type",
            ax=axes[1], show=False, title="RNA Velocity — Cell Types"
        )

    plt.tight_layout()
    plt.savefig("results/figures/velocity/velocity_streams.png", dpi=150)
    plt.close()
    log.info("  Saved: results/figures/velocity/velocity_streams.png")

    # Save
    out_path = "data/processed/GSE303787_velocity.h5ad"
    adata.write(out_path)
    log.info(f"\n✓ Velocity complete. Saved to: {out_path}")
    log.info("  Next: python src/scoring/engraftment_score.py")


if __name__ == "__main__":
    main()
