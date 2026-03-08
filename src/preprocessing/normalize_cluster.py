import scanpy as sc
import harmonypy as hm
import numpy as np
import matplotlib.pyplot as plt
import yaml
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

sc.settings.verbosity = 1
Path("results/figures/clustering").mkdir(parents=True, exist_ok=True)

CELL_TYPE_MARKERS = {
    "iPSC":    ["POU5F1", "NANOG", "SOX2", "LIN28A"],
    "NMP":     ["T", "SOX2", "CDX2", "MSGN1"],
    "spNPC":   ["SOX2", "NES", "VIM", "PAX6"],
    "pMN":     ["OLIG2", "NKX6-1", "DBX2"],
    "MN":      ["MNX1", "CHAT", "ISL1"],
    "V2_IN":   ["CHX10", "EN1", "FOXN4"],
    "Mesoderm":["MEOX1", "MEOX2", "TBX6"],
}

def main():
    cfg = yaml.safe_load(open("configs/qc_params.yaml"))

    # ── Step 1: Normalize ─────────────────────────────────────────────────────
    norm_path = "data/processed/GSE303787_normalized.h5ad"
    if Path(norm_path).exists():
        log.info(f"Loading checkpoint: {norm_path}")
        adata = sc.read_h5ad(norm_path)
    else:
        log.info("Loading QC-filtered data...")
        adata = sc.read_h5ad("data/processed/GSE303787_qc_filtered.h5ad")
        log.info(f"  {adata.n_obs} cells x {adata.n_vars} genes")

        log.info("Normalizing...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata

        log.info("Selecting HVGs...")
        sc.pp.highly_variable_genes(
            adata, n_top_genes=3000, flavor="seurat_v3",
            batch_key="sample_id", subset=False
        )
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata_hvg, max_value=10)

        log.info("Running PCA...")
        sc.tl.pca(adata_hvg, n_comps=50, svd_solver="arpack")

        adata_hvg.write(norm_path)
        log.info(f"  Checkpoint saved: {norm_path}")
        adata = adata_hvg

# ── Step 2: Harmony ───────────────────────────────────────────────────────
    harmony_path = "data/processed/GSE303787_harmony.h5ad"
    if Path(harmony_path).exists():
        log.info(f"Loading checkpoint: {harmony_path}")
        adata = sc.read_h5ad(harmony_path)
    else:
        log.info("Running Harmony batch correction...")
        ho = hm.run_harmony(
            adata.obsm["X_pca"], adata.obs, "sample_id",
            max_iter_harmony=30, random_state=42, verbose=False
        )
        Z = ho.Z_corr
        if Z.shape[0] == adata.n_obs:
            adata.obsm["X_pca_harmony"] = Z
        else:
            adata.obsm["X_pca_harmony"] = Z.T
        adata.write(harmony_path)
        log.info(f"  Checkpoint saved: {harmony_path}")

    # ── Step 3: Neighbors + Leiden ────────────────────────────────────────────
    leiden_path = "data/processed/GSE303787_leiden.h5ad"
    if Path(leiden_path).exists():
        log.info(f"Loading checkpoint: {leiden_path}")
        adata = sc.read_h5ad(leiden_path)
    else:
        log.info("Building neighborhood graph...")
        sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30,
                        use_rep="X_pca_harmony", random_state=42)

        log.info("Leiden clustering...")
        sc.tl.leiden(adata, resolution=0.5, random_state=42)
        n = adata.obs["leiden"].nunique()
        log.info(f"  {n} clusters found")

        adata.write(leiden_path)
        log.info(f"  Checkpoint saved: {leiden_path}")

    # ── Step 4: UMAP ──────────────────────────────────────────────────────────
    umap_path = "data/processed/GSE303787_clustered.h5ad"
    if Path(umap_path).exists():
        log.info(f"Loading checkpoint: {umap_path}")
        adata = sc.read_h5ad(umap_path)
    else:
        log.info("Computing UMAP (this is the slow step ~10-15 min)...")
        sc.tl.umap(adata, random_state=42)

        # Cell type scoring
        for ct, markers in CELL_TYPE_MARKERS.items():
            valid = [m for m in markers if m in adata.var_names]
            if valid:
                sc.tl.score_genes(adata, valid, score_name=f"score_{ct}",
                                  random_state=42)

        score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
        if score_cols:
            adata.obs["predicted_cell_type"] = (
                adata.obs[score_cols].idxmax(axis=1).str.replace("score_", "")
            )

        adata.write(umap_path)
        log.info(f"  Checkpoint saved: {umap_path}")

    # ── Step 5: Plot ──────────────────────────────────────────────────────────
    log.info("Plotting UMAP...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, title="Leiden Clusters")
    sc.pl.umap(adata, color="sample_id", ax=axes[1], show=False, title="Timepoint")
    if "predicted_cell_type" in adata.obs.columns:
        sc.pl.umap(adata, color="predicted_cell_type", ax=axes[2],
                   show=False, title="Cell Type")
    plt.tight_layout()
    plt.savefig("results/figures/clustering/umap_overview.png", dpi=150)
    plt.close()
    log.info("  Saved: results/figures/clustering/umap_overview.png")

    log.info("\n✓ Done. Next: python src/trajectory/rna_velocity.py")

if __name__ == "__main__":
    main()
