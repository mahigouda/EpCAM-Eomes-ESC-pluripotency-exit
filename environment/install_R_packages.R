## =============================================================================
## Install the R packages used by the scripts in this repository.
## Tested with R >= 4.2 and Bioconductor >= 3.16. Run once:
##   Rscript environment/install_R_packages.R
## =============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran_pkgs <- c(
  "tidyverse", "dplyr", "tidyr", "purrr", "readxl", "openxlsx", "WriteXLS", "gdata",
  "gtools", "knitr", "rmarkdown", "DT", "ggplot2", "ggrepel", "ggridges", "gplots",
  "pheatmap", "RColorBrewer", "viridis", "cowplot", "gridExtra", "gtable", "patchwork",
  "reshape2", "Matrix", "uwot", "VGAM", "princurve", "igraph", "Seurat"
)

bioc_pkgs <- c(
  "DESeq2", "IHW", "Rsubread", "pcaExplorer", "EnhancedVolcano", "AnnotationDbi",
  "EnsDb.Mmusculus.v79", "EnsDb.Hsapiens.v86", "org.Mm.eg.db", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot", "pathview", "SummarizedExperiment",
  "SingleCellExperiment", "scater", "scran", "slingshot", "TrajectoryUtils", "monocle",
  "Biobase", "BiocGenerics", "BiocParallel", "S4Vectors", "IRanges", "XVector",
  "GenomeInfoDb", "GenomicRanges", "GenomicAlignments", "Biostrings", "Rsamtools",
  "ShortRead", "MatrixGenerics", "matrixStats", "ExperimentHub", "depmap", "MAGeCKFlute"
)

github_pkgs <- c(
  "mojaveazure/seurat-disk",           # SeuratDisk
  "mojaveazure/loomR",                 # loomR
  "cole-trapnell-lab/monocle3"         # monocle3 (see its installation notes)
)

install.packages(setdiff(cran_pkgs, rownames(installed.packages())))
BiocManager::install(setdiff(bioc_pkgs, rownames(installed.packages())), update = FALSE, ask = FALSE)
for (repo in github_pkgs) {
  pkg <- switch(basename(repo), "seurat-disk" = "SeuratDisk", basename(repo))
  if (!requireNamespace(pkg, quietly = TRUE)) remotes::install_github(repo, upgrade = "never")
}

message("All packages installed. Record versions with sessionInfo().")
