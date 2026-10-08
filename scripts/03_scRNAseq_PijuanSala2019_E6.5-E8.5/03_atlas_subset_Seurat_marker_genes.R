## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 03 | Public scRNA-seq: mouse gastrulation atlas E6.5-E8.5
##           (Pijuan-Sala et al., Nature 2019; ArrayExpress E-MTAB-6967) and
##           mouse organogenesis E9.5-E11.5 (Dong et al., Genome Biol 2018; GEO: GSE87038)
## Script  : 03_atlas_subset_Seurat_marker_genes.R
## Purpose : Seurat SCTransform/PCA/UMAP of an atlas expression subset and feature plots of the
##           EpCAM-associated marker genes.
## Input   : data/PijuanSala2019_atlas/adata_2.txt (genes x cells table exported from the atlas AnnData)
## Output  : interactive plots
## =============================================================================

## ---- Packages -----------------------------------------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(scater)
library(monocle)
library(reshape2)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(viridis)
library(gridExtra)
library(ggridges)

## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/PijuanSala2019_atlas/ (see data/README.md); all outputs are written to
## results/03_PijuanSala2019/.
data_dir <- file.path(getwd(), "data", "PijuanSala2019_atlas")
out_dir  <- file.path(getwd(), "results", "03_PijuanSala2019")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- SCTransform, PCA (non-zero-variance genes), UMAP, marker genes ---------

data <- read.table(file.path(data_dir, "adata_2.txt"), header = TRUE, row.names = 1, sep = "\t")
sce <- CreateSeuratObject(counts = data)

# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)

# Calculate gene variance
gene_vars <- scater::perFeatureVariance(sce)

# Filter out genes with zero variance
non_zero_variance_genes <- rownames(gene_vars[gene_vars$variance != 0, ])

# Subset the SCE object to only include genes with non-zero variance
sce_filtered <- sce[non_zero_variance_genes, ]

# Run PCA on the filtered dataset
sce_filtered <- RunPCA(sce_filtered)


sce <- RunPCA(sce)

# Run UMAP
sce <- RunUMAP(sce, dims = 1:10)

# Visualize the results
DimPlot(sce, reduction = "umap")

genes_of_interest <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112", "Nanog-6-122707489",
                       "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
                       "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Fgf4-7-144861386", "Esrrb-12-86361117",
                       "Twist1-12-33957671", "Cdx2-5-147300805", "Hand1-11-57828705", "T-17-8434423" )

FeaturePlot(sce, features = genes_of_interest, reduction = "umap")
