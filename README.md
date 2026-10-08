# EpCAM and exit from pluripotency: transcriptomic analyses

![R](https://img.shields.io/badge/R-DESeq2%20%7C%20Seurat%20%7C%20Monocle3-276DC3?logo=r)
![Python](https://img.shields.io/badge/Python-Scanpy-3776AB?logo=python)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Analysis scripts for bulk and single-cell RNA-seq data on the role of the surface protein
EpCAM during the differentiation of mouse embryonic stem cells (mESC). The scripts cover
bulk RNA-seq of wild-type and *Epcam*-knockout mESC and the re-analysis of public
single-cell RNA-seq atlases of mouse, human and zebrafish embryos.

## Related publication

The biological context of these analyses is described in:

> Gong N, Gouda M, *et al.* EpCAM supports exit from pluripotency of embryonic stem cells via
> Eomes. *Cell Death & Disease* 17, 389 (2026). <https://doi.org/10.1038/s41419-026-08734-w>

Please refer to the article for the study design, experimental methods and the validated
results.

## Overview of analyses

| Module | Data | Methods |
|---|---|---|
| `01_bulk_RNAseq` | 3′-RNA-seq of wild-type and two *Epcam*<sup>−/−</sup> mESC clones during embryoid-body differentiation (D0–D10), GEO GSE293121 | DESeq2 with IHW, PCA, volcano plots, heatmaps, GO enrichment (clusterProfiler), candidate-gene expression |
| `02_scRNAseq_Mohammed2017_E3.5-E6.75` | Mouse peri-implantation embryos (Mohammed *et al.* 2017), GEO GSE100597 | Seurat, lineage module scores, gene co-expression and correlation, Monocle3 and Slingshot trajectories |
| `03_scRNAseq_PijuanSala2019_E6.5-E8.5` | Mouse gastrulation atlas (Pijuan-Sala *et al.* 2019), E-MTAB-6967; mouse organogenesis (Dong *et al.* 2018), GEO GSE87038 | Scanpy preprocessing, stage-level pseudobulk, Seurat analysis of heart cells |
| `04_scRNAseq_human_embryo_CS12-CS16` | Human embryos (Xu *et al.* 2023), GEO GSE157329 | Scanpy, marker expression by stage and organ system, single-cell co-expression |
| `05_scRNAseq_zebrafish` | Zebrafish lineage atlas (Farrell *et al.* 2018), GEO GSE106587 | Scanpy, lineage annotation, marker expression and co-expression |

## Repository structure

```
scripts/
├── 01_bulk_RNAseq/
├── 02_scRNAseq_Mohammed2017_E3.5-E6.75/
├── 03_scRNAseq_PijuanSala2019_E6.5-E8.5/
├── 04_scRNAseq_human_embryo_CS12-CS16/
└── 05_scRNAseq_zebrafish/
data/README.md     # data sources and accession numbers
environment/       # R package list and conda environment
```

Each script begins with a short header that describes its purpose, inputs and outputs.
The Python scripts use the jupytext "percent" format and can be opened as notebooks.

## Data

No data are included in this repository. All datasets are publicly available. Sources and
accession numbers are listed in [`data/README.md`](data/README.md).

## Software

R (Bioconductor: DESeq2, IHW, clusterProfiler, SingleCellExperiment, slingshot; CRAN:
Seurat, tidyverse, pheatmap) and Python (Scanpy, pandas, NumPy, SciPy, matplotlib, seaborn).
See [`environment/`](environment/).

## Disclaimer

These scripts are shared for reference and transparency. They are **not** an official or
maintained software release of the publication above, and they are not a packaged,
end-to-end pipeline. They were developed interactively in an exploratory research setting
and are provided as is. Results may differ from the published figures because of differences
in software versions, random seeds, data versions or preprocessing. For the authoritative
methods and results, please refer to the published article.

## Citation

If these scripts are useful for your work, please cite the related article (see above and
[`CITATION.cff`](CITATION.cff)) and the original publications of the public datasets listed in
[`data/README.md`](data/README.md).

## License

Code is available under the [MIT License](LICENSE). Datasets remain subject to the terms of
their original repositories.
