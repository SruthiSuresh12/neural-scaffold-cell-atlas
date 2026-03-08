# neural-scaffold-cell-atlas

**Computational characterization of iPSC-derived spinal cord progenitor cell states to define optimal engraftment criteria for biomaterial scaffolds**

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Motivation

Scaffold-mediated cell therapy for spinal cord injury depends critically on transplanting iPSC-derived progenitors that will survive and integrate into host tissue. Current protocols select cells based on bulk markers alone — ignoring transcriptional heterogeneity that determines whether individual cells will engraft successfully or undergo apoptosis post-transplant.

This project builds a single-cell transcriptomic atlas of human iPSC-derived spinal cord progenitors, reconstructs developmental trajectories, and produces a per-cell **Engraftment Suitability Score (ESS)** — a molecular readout predicting which cell states are most likely to survive transplantation into a biomaterial scaffold.

## Pipeline
```
Raw scRNA-seq (GEO)
      │
      ▼
1. Quality Control        ← doublet removal, MT% filtering
      │
      ▼
2. Normalization          ← HVG selection, PCA, Harmony batch correction
      │
      ▼
3. Clustering             ← Leiden, UMAP, cell type annotation
      │
      ▼
4. Trajectory Inference   ← scVelo RNA velocity
      │
      ▼
5. Engraftment Score      ← ESS: survival markers vs. stress signatures
```

## Dataset

| GEO Accession | Description | Reference |
|---|---|---|
| [GSE303787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE303787) | Human iPSC → spinal cord progenitors, Day 1–18 | Yao, Zhang et al., *Cell Death & Disease* 2025 ([PMID: 40846836](https://pubmed.ncbi.nlm.nih.gov/40846836/)) |
| [GSE171890](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171890) | Human fetal spinal cord reference (CS12–CS19) | Rayon et al., *Development* 2021 ([DOI: 10.1242/dev.199711](https://doi.org/10.1242/dev.199711)) |

## Key Results

**Optimal harvest window:** Day18 spNeuron population achieves highest median ESS (0.626) with 37.5% of cells in the high-suitability tier — consistent with maximal neural maturity markers at this timepoint.

**ESS by timepoint:**
| Timepoint | Median ESS | % High-Tier |
|---|---|---|
| Day18_spNeuron | 0.626 | 37.5% |
| Day6_NMP | 0.510 | 9.1% |
| Day12_spNPG | 0.502 | 11.2% |
| Day3_NMP | 0.422 | 0.6% |
| Day1_iPSC | 0.400 | 0.4% |


## Installation
```bash
git clone https://github.com/SruthiSuresh12/neural-scaffold-cell-atlas.git
cd neural-scaffold-cell-atlas
conda env create -f environment.yml
conda activate scaffold-atlas
```

## Usage
```bash
# Download data from GEO GSE303787, organize into data/raw/GSE303787/{Day1,Day3,Day6,Day12,Day18}/
# Then run pipeline step by step:
python src/preprocessing/qc_filter.py
python src/preprocessing/normalize_cluster.py
python src/trajectory/rna_velocity.py
python src/scoring/engraftment_score.py

# Or end-to-end:
snakemake --cores 4
```

## References

- Yao X et al. (2025). Human iPSC-derived spinal neural progenitors enhance sensorimotor recovery in spinal cord-injured NOD-SCID mice. *Cell Death & Disease*. [DOI: 10.1038/s41419-025-07961-x](https://doi.org/10.1038/s41419-025-07961-x)
- Chew SY et al. (2025). Sustained MicroRNA delivery enhanced remyelination after spinal cord injury. *Biomaterials*. [DOI: 10.1016/j.biomaterials.2025.123565](https://doi.org/10.1016/j.biomaterials.2025.123565)
- Chew SY et al. (2025). Bioprinted Microchannel Scaffolds Modulate Neuronal Differentiation. *ACS Applied Bio Materials*. [DOI: 10.1021/acsabm.5c00441](https://doi.org/10.1021/acsabm.5c00441)
- Rayon T et al. (2021). Single-cell transcriptome profiling of the human developing spinal cord. *Development*. [DOI: 10.1242/dev.199711](https://doi.org/10.1242/dev.199711)

## Author

**Sruthi Suresh** | B.Tech CSE (Bioinformatics), VIT | Visiting Research Intern, EPFL Brain Mind Institute
[GitHub](https://github.com/SruthiSuresh12) | [LinkedIn](https://linkedin.com/in/sruthi-suresh-sruthi-suresh)
EOF
