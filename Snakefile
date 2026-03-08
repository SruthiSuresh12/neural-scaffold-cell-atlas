rule all:
    input:
        "data/processed/GSE303787_ESS_scored.h5ad",
        "results/figures/scoring/ESS_panel_figure.png",
        "results/tables/optimal_harvest_timepoints.csv"

rule qc_filter:
    output: "data/processed/GSE303787_qc_filtered.h5ad"
    shell: "python src/preprocessing/qc_filter.py"

rule normalize_cluster:
    input: "data/processed/GSE303787_qc_filtered.h5ad"
    output: "data/processed/GSE303787_clustered.h5ad"
    shell: "python src/preprocessing/normalize_cluster.py"

rule rna_velocity:
    input: "data/processed/GSE303787_clustered.h5ad"
    output: "data/processed/GSE303787_velocity.h5ad"
    shell: "python src/trajectory/rna_velocity.py"

rule engraftment_score:
    input: "data/processed/GSE303787_clustered.h5ad"
    output:
        "data/processed/GSE303787_ESS_scored.h5ad",
        "results/figures/scoring/ESS_panel_figure.png",
        "results/tables/optimal_harvest_timepoints.csv"
    shell: "python src/scoring/engraftment_score.py"
