# Transcriptomic analysis of EpCAM-dependent exit from pluripotency

**Bulk RNA-seq · single-cell RNA-seq · cross-species embryo atlases · R & Python**

[![Publication](https://img.shields.io/badge/Published%20in-Cell%20Death%20%26%20Disease%20(2026)-blue)](https://doi.org/10.1038/s41419-026-08734-w)
![R](https://img.shields.io/badge/R-DESeq2%20%7C%20Seurat%20%7C%20Monocle3-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-Scanpy%20%7C%20pandas-3776AB?logo=python)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

This repository is a portfolio of the bioinformatic analyses I carried out in a stem-cell
research project on how the surface protein **EpCAM** helps mouse embryonic stem cells (mESC)
leave the pluripotent state and differentiate. The biological results were published in
*Cell Death & Disease* (2026), and I am a co-author of the article:

> Gong N, Gouda M, *et al.* EpCAM supports exit from pluripotency of embryonic stem cells via
> Eomes. *Cell Death & Disease* 17, 389 (2026). <https://doi.org/10.1038/s41419-026-08734-w>

This repository is not the official code release of the article. It shows how I approached
the analyses end to end: designing differential-expression contrasts, integrating public
single-cell atlases from three species, and turning results into publication figures.

---

## The question

In stem cells that lack EpCAM, differentiation into heart muscle (cardiomyocytes) fails.
I wanted to find out which genes and pathways explain this, and whether the candidate genes
are also active together with EpCAM in real embryos.

## What I did

| Step | Data | What I analysed | Key tools |
|---|---|---|---|
| **1. Bulk RNA-seq: wild type vs knockout** | WT and two *Epcam*-knockout clones × 4 differentiation days (D0–D10) × 4 replicates | 21 differential-expression contrasts (clone vs WT, pooled knockout vs WT, time courses), PCA, volcano plots, heatmaps, GO enrichment, candidate genes (Eomes, Foxa2, Gata6) | DESeq2, IHW, EnhancedVolcano, clusterProfiler, pheatmap |
| **2. Early mouse embryo (E3.5–E6.75)** | Public Smart-seq2 scRNA-seq (Mohammed *et al.* 2017) | Cell-lineage assignment from module scores, *Epcam*/*Eomes* co-expression, per-stage gene–gene correlation, pseudotime trajectories | Seurat (SCTransform), Monocle3, Slingshot |
| **3. Mouse gastrulation & organogenesis (E6.5–E11.5)** | Public atlases: Pijuan-Sala *et al.* 2019 (~116k cells), Dong *et al.* 2018 | Large-scale preprocessing, stage-level pseudobulk, heart-cell lineage scoring | Scanpy, Seurat |
| **4. Human embryo (Carnegie stages 12–16)** | Public 10x scRNA-seq (Xu *et al.* 2023, ~185k cells) | Marker expression by stage and organ system, single-cell co-expression of gene combinations | Scanpy, seaborn |
| **5. Zebrafish embryo** | Public Drop-seq lineage atlas (Farrell *et al.* 2018) | Lineage annotation, marker expression and co-expression across stages | Scanpy |

**Main finding supported by these analyses:** loss of EpCAM changes Wnt-signalling and
heart-development genes early in differentiation. Of the candidates, *Eomes* stood out, and
it is co-expressed with *Epcam* in early-lineage cells of mouse and human embryos.

## Skills shown

- **Experimental design for RNA-seq:** multi-factor contrasts (genotype × time), replicate
  handling, testing against a fold-change threshold with independent hypothesis weighting
- **Single-cell analysis:** QC, normalisation (SCTransform / log-normalisation), dimensionality
  reduction, clustering, module-score cell annotation, trajectory inference
- **Data integration:** re-using public data across mouse, human and zebrafish, and atlases of
  100k+ cells
- **Two languages:** R / Bioconductor and Python / Scanpy in one project
- **Communication:** publication-quality figures (volcano plots, dot plots, UMAPs, heatmaps,
  correlation matrices) for a peer-reviewed article

## Repository structure

```
scripts/
├── 01_bulk_RNAseq/                        # DESeq2 report, DE contrasts, volcano plots, GO
├── 02_scRNAseq_Mohammed2017_E3.5-E6.75/   # Seurat lineage scoring, co-expression, trajectories
├── 03_scRNAseq_PijuanSala2019_E6.5-E8.5/  # Scanpy atlas preprocessing, heart-cell analysis
├── 04_scRNAseq_human_embryo_CS12-CS16/    # human embryo co-expression (Scanpy)
└── 05_scRNAseq_zebrafish/                 # zebrafish lineage atlas (Scanpy)
data/README.md                             # public data sources and accession numbers
environment/                               # R package list and conda environment
```

Each script starts with a short header that explains its purpose, inputs and outputs. The
Python scripts use the jupytext "percent" format, so they open as notebooks in Jupyter or VS Code.

## Data

No data are included. All datasets are public. Accession numbers and download links are in
[`data/README.md`](data/README.md): GEO GSE293121, GSE100597, GSE87038, GSE157329 and
GSE106587, and ArrayExpress E-MTAB-6967.

## Note on the code

The scripts come from an exploratory research project and were developed interactively in
RStudio and Jupyter. I organised them by analysis and added headers, but kept the analysis
code as I wrote it. They are shared to show methods and workflow, not as a packaged pipeline.

## About me

I am a bioinformatician working on transcriptomics and single-cell genomics.
Feel free to connect on [LinkedIn](https://www.linkedin.com/in/<your-profile>) or open an
issue in this repository.

## License

Code: [MIT](LICENSE). Datasets remain under the terms of their original repositories.
