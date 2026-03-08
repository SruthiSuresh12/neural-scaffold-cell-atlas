import scanpy as sc
import scrublet as scr
import numpy as np
import matplotlib.pyplot as plt
import yaml
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

sc.settings.verbosity = 1
Path("results/figures/qc").mkdir(parents=True, exist_ok=True)
Path("data/processed").mkdir(parents=True, exist_ok=True)


def load_config():
    with open("configs/qc_params.yaml") as f:
        return yaml.safe_load(f)


def load_sample(sample_dir, sample_name):
    import scipy.io
    import pandas as pd
    log.info(f"Loading {sample_name} from {sample_dir}")

    matrix = scipy.io.mmread(str(sample_dir / "matrix.mtx.gz")).T.tocsr()

    barcodes = pd.read_csv(
        str(sample_dir / "barcodes.tsv.gz"), header=None
    )[0].values

    feat_path = sample_dir / "features.tsv.gz"
    genes_df = pd.read_csv(str(feat_path), header=None, sep="\t")
    gene_names = genes_df[1].values

    import anndata
    adata = anndata.AnnData(
        X=matrix,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=gene_names)
    )
    adata.obs_names = [f"{sample_name}_{bc}" for bc in adata.obs_names]
    adata.obs["sample_id"] = sample_name
    adata.var_names_make_unique()
    log.info(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    return adata

def qc_metrics(adata):
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["rb"] = adata.var_names.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rb"], inplace=True)
    return adata


def doublet_detection(adata):
    log.info("  Running Scrublet...")
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.06)
    scores, predicted = scrub.scrub_doublets(verbose=False)
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted
    log.info(f"  Doublets detected: {predicted.sum()} ({100*predicted.sum()/adata.n_obs:.1f}%)")
    return adata


def filter_cells(adata, cfg):
    qc = cfg["qc"]
    n_before = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    adata = adata[
        (adata.obs["n_genes_by_counts"] >= qc["min_genes"]) &
        (adata.obs["n_genes_by_counts"] <= qc["max_genes"]) &
        (adata.obs["pct_counts_mt"] <= qc["max_mt_percent"]) &
        (adata.obs["pct_counts_rb"] <= qc["max_rb_percent"])
    ].copy()
    sc.pp.filter_genes(adata, min_cells=qc["min_cells"])
    log.info(f"  Cells: {n_before} -> {adata.n_obs} retained")
    return adata


def main():
    cfg = load_config()

    samples = {
        "Day1_iPSC":      Path("data/raw/GSE303787/Day1"),
        "Day3_NMP":       Path("data/raw/GSE303787/Day3"),
        "Day6_NMP":       Path("data/raw/GSE303787/Day6"),
        "Day12_spNPG":    Path("data/raw/GSE303787/Day12"),
        "Day18_spNeuron": Path("data/raw/GSE303787/Day18"),
    }

    adatas = []
    for name, path in samples.items():
        if not path.exists():
            log.warning(f"  Skipping {name} — path not found: {path}")
            continue
        adata = load_sample(path, name)
        adata = qc_metrics(adata)
        adata = doublet_detection(adata)
        adata = filter_cells(adata, cfg)
        adatas.append(adata)

    if not adatas:
        log.error("No samples loaded. Check data paths.")
        return

    log.info(f"\nConcatenating {len(adatas)} samples...")
    combined = adatas[0].concatenate(
        adatas[1:],
        batch_key="sample_id",
        batch_categories=list(samples.keys())[:len(adatas)]
    )
    log.info(f"Combined: {combined.n_obs} cells x {combined.n_vars} genes")

    out = "data/processed/GSE303787_qc_filtered.h5ad"
    combined.write(out)
    log.info(f"\n✓ Saved to {out}")
    log.info("Next step: python src/preprocessing/normalize_cluster.py")


if __name__ == "__main__":
    main()

