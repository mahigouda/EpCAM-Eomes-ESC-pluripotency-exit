## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 02 | Public scRNA-seq: mouse peri-implantation embryos E3.5-E6.75
##           (Mohammed et al., Cell Rep 2017; GEO: GSE100597)
## Script  : 04_Slingshot_trajectory.R
## Purpose : Lineage/pseudotime inference with Slingshot on the E3.5-E6.75 embryo cells.
## Input   : data/Mohammed2017_GSE100597/GSE100597_count_table_QC_filtered.txt
## Output  : results/02_Mohammed2017/GSE100597_count.rds (+ interactive plots)
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
library(loomR)
library(princurve)
library(TrajectoryUtils)
library(slingshot)
library(ggrepel)

## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/Mohammed2017_GSE100597/ (see data/README.md); all outputs are written to
## results/02_Mohammed2017/.
data_dir <- file.path(getwd(), "data", "Mohammed2017_GSE100597")
out_dir  <- file.path(getwd(), "results", "02_Mohammed2017")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Load counts, SCTransform, UMAP, clustering ----------------------------

# Load the data
data <- read.table(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header = TRUE, row.names = 1, sep = "\t")
sce <- CreateSeuratObject(counts = data)

# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)

# Prepare data for Slingshot
sce <- RunUMAP(sce, dims = 1:10) # if you haven't run UMAP yet
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5) # adjust resolution as necessary


## ---- Slingshot on Seurat UMAP/clusters ------------------------------------
## Pseudotime analysis
# Extract cluster labels and reduced dimensional coordinates from Seurat object
rd <- sce@reductions$umap@cell.embeddings
cl <- sce@meta.data$seurat_clusters
# Prepare a SingleCellExperiment object
sce_sling <- SingleCellExperiment(list(counts = sce@assays$RNA@counts), colData = data.frame(cl))
# Add reduced dimension data to the SingleCellExperiment object
reducedDims(sce_sling) <- SimpleList(UMAP = rd)
# Set row data
rowData(sce_sling)$feature_symbol <- rownames(sce_sling)
# Run Slingshot
sds <- slingshot(sce_sling, clusterLabels = 'cl', reducedDim = 'UMAP')
# Plot the results
plot(reducedDims(sds)$UMAP, col = brewer.pal(9,'Set1')[sce@meta.data$seurat_clusters])
lines(SlingshotDataSet(sds), lwd = 2, type = 'lineages')

## Improve Pseudotime plots
# Create a data frame from UMAP coordinates and cluster labels
df <- data.frame(UMAP1 = reducedDims(sds)$UMAP[,1], UMAP2 = reducedDims(sds)$UMAP[,2],
                 cluster = sce@meta.data$seurat_clusters)
# Add cluster centroids
centroids <- aggregate(cbind(UMAP1, UMAP2) ~ cluster, data = df, FUN = function(x) mean(x))

# Get pseudotime values
pseudotime <- slingPseudotime(sds)

# Add to dataframe
df$pseudotime <- pseudotime[,1]  # if you have multiple lineages, adjust the index accordingly

# Create the trajectory data
trajectory_data <- df[order(df$pseudotime),]

# Plot the results with ggplot2
ggplot(df, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
  geom_point(size = 1) +
  geom_path(data = trajectory_data, aes(group = cluster), color = 'black') +
  scale_color_brewer(palette = "Set1") +
  geom_text_repel(data = centroids, aes(x = UMAP1, y = UMAP2, label = cluster)) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(color = "Cluster")


# Extract pseudotime values for each lineage
pseudotime_values <- slingPseudotime(sds)

# Add them to the dataframe
for (i in 1:ncol(pseudotime_values)) {
  df[, paste0("pseudotime_", i)] <- pseudotime_values[, i]
}

# Create a list to store the plots
plot_list <- list()

# Loop through each lineage
for (i in 1:ncol(pseudotime_values)) {
  # Subset the data for each lineage
  lineage_data <- df[!is.na(df[, paste0("pseudotime_", i)]), ]

  # Order by pseudotime
  lineage_data <- lineage_data[order(lineage_data[, paste0("pseudotime_", i)]), ]

  # Create a plot for each lineage
  plot_list[[i]] <- ggplot(lineage_data, aes(x = UMAP1, y = UMAP2, color = factor(cluster))) +
    geom_point(size = 1) +
    geom_path(aes(group = cluster), color = 'black') +
    scale_color_brewer(palette = "Set1") +
    geom_text_repel(data = centroids, aes(x = UMAP1, y = UMAP2, label = cluster)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(color = "Cluster", title = paste("Lineage", i))
}

# Print the plots
for (i in 1:length(plot_list)) {
  print(plot_list[[i]])
}


## ---- Slingshot on LogNormalize workflow (getLineages/getCurves) ----------
## Slingshot
data <- read.delim(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header = T, row.names = 1)
# Extract characters before the first underscore in each column name
colnames(data) <- sub("_.*", "", colnames(data))
# Print the modified data with updated column names
print(data)

comp_matrix <- Matrix::Matrix(as.matrix(data), sparse = T)
saveRDS(comp_matrix, "GSE100597_count.rds")

umi_counts <- readRDS("GSE100597_count.rds")
dim(umi_counts)

# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

suppressPackageStartupMessages({
  library(scran)
  library(igraph)
})

suppressPackageStartupMessages({
})

# Data analysis with Seurat pipeline
data <- CreateSeuratObject(counts = umi_counts)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)

data <- FindNeighbors(data)
data <- FindClusters(data, resolution = 1)

data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)

# Plot the clusters
DimPlot(data, group.by = "RNA_snn_res.1")

# Save the objects as separate matrices for input in slingshot
dimred <- data@reductions$umap@cell.embeddings
clustering <- data$RNA_snn_res.1
counts <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])

suppressPackageStartupMessages({
})

# Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)

lineages

# Plot the lineages
par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
  text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(lineages, lwd = 3, col = "black")

#Run default Slingshot
set.seed(1)
lineages <- getLineages(data = dimred,
                        clusterLabels = clustering,
                        #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                        start.clus = "0") #define where to start the trajectories

lineages

#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){
  text( mean(dimred[clustering==i,1]),
        mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],  pch = 16)
lines(lineages, lwd = 3, col = 'black')

curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(curves, lwd = 3, col = "black")
