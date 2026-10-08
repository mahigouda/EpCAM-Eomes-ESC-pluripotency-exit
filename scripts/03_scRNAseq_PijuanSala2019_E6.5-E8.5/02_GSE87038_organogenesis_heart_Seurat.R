## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 03 | Public scRNA-seq: mouse gastrulation atlas E6.5-E8.5
##           (Pijuan-Sala et al., Nature 2019; ArrayExpress E-MTAB-6967) and
##           mouse organogenesis E9.5-E11.5 (Dong et al., Genome Biol 2018; GEO: GSE87038)
## Script  : 02_GSE87038_organogenesis_heart_Seurat.R
## Purpose : Seurat analysis of mouse organogenesis scRNA-seq (E9.5-E11.5 tissues) with focus on
##           the heart: expression of EpCAM-associated genes, lineage module scores (incl.
##           cardiomyocytes) and Epcam/Gata6/Wnt11 expression on lineage UMAPs.
## Input   : data/Dong2018_GSE87038/{GSE87038_Mouse_Organogenesis_UMI_counts_matrix.txt, heart.txt}
## Output  : results/03_GSE87038_organogenesis/ (normalised expression tables, gene expression csv)
## Notes   : heart.txt = heart-derived cells of GSE87038 (genes x cells, first column = gene).
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
library(SingleCellExperiment)
library(uwot)
library(Matrix)


## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/Dong2018_GSE87038/ (see data/README.md); all outputs are written to
## results/03_GSE87038_organogenesis/.
data_dir <- file.path(getwd(), "data", "Dong2018_GSE87038")
out_dir  <- file.path(getwd(), "results", "03_GSE87038_organogenesis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Whole dataset: QC, LogNormalize, PCA, marker genes -------------------
data <- read.table(file.path(data_dir, "GSE87038_Mouse_Organogenesis_UMI_counts_matrix.txt"), header=T, row.names=1, sep="\t")
head(data)

df <-  CreateSeuratObject(counts = data, project = "GSE87038", min.cells = 3, min.features = 200)

# Visualize QC metrics as a violin plot
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Increase the time limit for subset() function
options(expressions = 10000)

# Subset cells based on quality control criteria
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- NormalizeData(df)


# Save normalized expression data as a text file
write.table(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expression_GSE87038_Mouse_Organogenesis.txt", sep = "\t", quote = FALSE)

# Save normalized expression data as a CSV file
write.csv(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expressionGSE87038_Mouse_Organogenesis.csv", row.names = TRUE)


df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling the data
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)

## Perform linear dimensional reduction
df <- RunPCA(df, features = VariableFeatures(object = df))

# Examine and visualize PCA results a few different ways
print(df[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(df, dims = 1:2, reduction = "pca")

DimPlot(df, reduction = "pca")

## Overall marker genes
genes_of_interest <- c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog",
                       "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
                       "Stat3", "Fgf8", "Fgf10", "Fgf4", "Esrrb",
                       "Twist1", "Cdx2", "Hand1", "T" )

FeaturePlot(df, features = genes_of_interest, reduction = "pca")


expression_data <- as.data.frame(GetAssayData(object = df, slot = "data")[genes_of_interest, ])

write.csv(expression_data, file = "analysis_expressionGSE87038_Mouse_Organogenesis.csv")


# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature
# in each cluster. The color represents the average expression level
DotPlot(df, features = genes_of_interest) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(df, downsample = 100), features = genes_of_interest, size = 3)


## ---- Heart subset (heart.txt): LogNormalize, PCA/UMAP, module score -------
## Heart cluster
data <- read.table(file.path(data_dir, "heart.txt"), header=T, sep="\t")

df <-  CreateSeuratObject(counts = data, project = "GSE87038", min.cells = 3, min.features = 200)

# Visualize QC metrics as a violin plot
VlnPlot(df, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")


# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(df, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Increase the time limit for subset() function
options(expressions = 10000)

# Subset cells based on quality control criteria
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)

df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- NormalizeData(df)


# Save normalized expression data as a text file
write.table(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expression_heart_Mouse_Organogenesis.txt", sep = "\t", quote = FALSE)

# Save normalized expression data as a CSV file
write.csv(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expression_heart_Mouse_Organogenesis.csv", row.names = TRUE)


df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(df), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(df)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling the data
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)

## Perform linear dimensional reduction
df <- RunPCA(df, features = VariableFeatures(object = df))

# Examine and visualize PCA results a few different ways
print(df[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(df, dims = 1:2, reduction = "pca")

DimPlot(df, reduction = "pca")


# Compute UMAP embedding
object <- RunUMAP(df, dims = 1:20, metric = "cosine")

# Visualize UMAP projection
DimPlot(object, reduction = "umap", label = TRUE, pt.size = 1)

# Visualize expression of multiple genes
FeaturePlot(object, features = c("Epcam", "Foxa2", "Hand1"), pt.size = 1)


## Overall marker genes
genes_of_interest <- c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog",
                       "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
                       "Stat3", "Fgf8", "Fgf10", "Fgf4", "Esrrb",
                       "Twist1", "Cdx2", "Hand1", "T" )

# Set the minimum expression cutoff for all genes in genes_of_interest to 0
object <- AddModuleScore(object = df,
                         features = genes_of_interest,
                         set.ident = TRUE,
                         min.cells = 0)


FeaturePlot(object, features = genes_of_interest, reduction = "pca")

genes_of_interest <- c("Fgf8", "Fgf10", "Fgf4")

FeaturePlot(df, features = genes_of_interest, reduction = "pca")


expression_data <- as.data.frame(GetAssayData(object = df, slot = "data")[genes_of_interest, ])

write.csv(expression_data, file = "analysis_expression_heart_Mouse_Organogenesis.csv")


## ---- Whole dataset: SCTransform, UMAP, per-cluster sub-UMAPs ---------------
## New analysis
data <- read.table(file.path(data_dir, "GSE87038_Mouse_Organogenesis_UMI_counts_matrix.txt"), header=T, row.names=1, sep="\t")

sce <- CreateSeuratObject(counts = data)

# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)
# Run UMAP
sce <- RunUMAP(sce, dims = 1:10)

# Visualize the results
DimPlot(sce, reduction = "umap")

genes_of_interest <- c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog",
                       "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
                       "Stat3", "Fgf8", "Fgf10", "Fgf4", "Esrrb",
                       "Twist1", "Cdx2", "Hand1", "T" )

FeaturePlot(sce, features = genes_of_interest, reduction = "umap")

# 1. Identify clusters
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.5)  # Adjust the resolution parameter according to your dataset

# 2. Subset the main Seurat object based on the identified clusters
clusters <- unique(sce$seurat_clusters)
subsets <- lapply(clusters, function(cluster) {
  return(subset(sce, subset = seurat_clusters == cluster))
})

# 3. Run UMAP on each subset and store the results in a list
sub_umaps <- lapply(subsets, function(subset_sce) {
  subset_sce <- NormalizeData(subset_sce)
  subset_sce <- ScaleData(subset_sce)
  subset_sce <- RunPCA(subset_sce)
  subset_sce <- RunUMAP(subset_sce, dims = 1:10)
  return(subset_sce)
})

# Optional: Combine the UMAP results into a single plot

plots <- lapply(sub_umaps, function(subset_sce) {
  plot <- DimPlot(subset_sce, group.by = "seurat_clusters", label = TRUE, label.size = 4)
  return(plot)
})

combined_plot <- wrap_plots(plots, ncol = 2)  # Adjust ncol according to the number of clusters
combined_plot


head(sce@meta.data)

colnames(sce@meta.data)
rownames(sce@meta.data)


## ---- Tissue subsets (heart, liver) -------------------------------------
## Heart
heart_indices <- grep("heart", rownames(sce@meta.data))
sce_heart <- sce[heart_indices, ]

sce_heart <- NormalizeData(sce_heart)

sce_heart <- FindVariableFeatures(sce_heart)

sce_heart <- ScaleData(sce_heart)
sce_heart <- RunPCA(sce_heart)

sce_heart <- RunUMAP(sce_heart, dims = 1:10)

sce_heart <- FindNeighbors(sce_heart, dims = 1:10)
sce_heart <- FindClusters(sce_heart, resolution = 0.5)

DimPlot(sce_heart, group.by = "seurat_clusters")

head(sce_heart@meta.data)

## liver
liver_indices <- grep("liver", rownames(sce@meta.data))
sce_liver <- sce[liver_indices, ]

sce_liver <- NormalizeData(sce_liver)

sce_liver <- FindVariableFeatures(sce_liver)

sce_liver <- ScaleData(sce_liver)
sce_liver <- RunPCA(sce_liver)

sce_liver <- RunUMAP(sce_liver, dims = 1:10)

sce_liver <- FindNeighbors(sce_liver, dims = 1:10)
sce_liver <- FindClusters(sce_liver, resolution = 0.5)

DimPlot(sce_liver, group.by = "seurat_clusters")

head(sce_liver@meta.data)

DimPlot(sce_liver, group.by = "SCT_snn_res.0.5")


liver_indices <- which(sce@meta.data$orig.ident == "liver")
sce_liver <- sce[liver_indices, ]
sce_liver <- RunPCA(sce_liver)
sce_liver <- RunUMAP(sce_liver, dims = 1:10)
DimPlot(sce_liver, group.by = "orig.ident")


## ---- Heart subset: SCTransform, lineage module scores, Epcam/Gata6/Wnt11 ----
data <- read.table(file.path(data_dir, "heart.txt"), header = TRUE, sep = "\t")
rownames(data) <- make.unique(data[, 1])
data <- data[, -1]

sce <- CreateSeuratObject(counts = data)
# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)

# Run UMAP
sce <- RunUMAP(sce, dims = 1:10)
DimPlot(sce, reduction = "umap")

## Marker genes expression
genes_of_interest <- c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog",
                       "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
                       "Fgf8", "Fgf10", "Fgf4", "Esrrb",
                       "Twist1", "Cdx2", "Hand1", "T" )

FeaturePlot(sce, features = genes_of_interest, reduction = "umap")

## Cell type assigning and module scores
known_markers <- list(
  pluripotency = c("Sox2", "Nanog", "Lin28a"),
  mesoendoderm = c("Eomes"),
  mesoderm = c("Tbx6","Hand1"),
  endoderm = c("Sox17", "Sox7", "Gata4"),
  ectoderm = c("Pax6", "Neurod1", "Nes"),
  cardiomyocytes = c("Myl4", "Tnnt2", "Myl7", "Tnni1", "Nexn",
                     "Acta2", "Mef2c", "Actc1", "Gyg"),
  epiblast = c("Dnmt3b","Gng3",
               "Gstt2", "Dnmt3a"))


# Calculate the module scores for each cell type
for (cell_type in names(known_markers)) {
  sce <- AddModuleScore(sce, features = known_markers[[cell_type]], name = cell_type)
}

module_score_names <- paste0("ModuleScore_", names(known_markers))

# Create a list of module score names for each cell type
module_score_names_list <- lapply(names(known_markers), function(cell_type) {
  prefix <- paste0(cell_type, "\\d+")
  matching_cols <- grep(prefix, colnames(sce@meta.data), perl = TRUE, value = TRUE)
  return(matching_cols)
})

# Calculate the average module score for pluripotency
pluripotency_scores <- sce@meta.data[, grep("pluripotency\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_pluripotency <- rowMeans(pluripotency_scores)

# Calculate the module score for mesoendoderm
mesoendoderm_scores <- sce@meta.data[, grep("mesoendoderm\\d+", colnames(sce@meta.data), perl = TRUE)]
# Add the individual module scores to the meta.data object
sce@meta.data <- cbind(sce@meta.data, mesoendoderm_scores)


# Calculate the average module score for mesoderm
mesoderm_scores <- sce@meta.data[, grep("mesoderm\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_mesoderm <- rowMeans(mesoderm_scores)

# Calculate the average module score for endoderm
endoderm_scores <- sce@meta.data[, grep("endoderm\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_endoderm <- rowMeans(endoderm_scores)

# Calculate the average module score for ectoderm
ectoderm_scores <- sce@meta.data[, grep("ectoderm\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_ectoderm <- rowMeans(ectoderm_scores)

# Calculate the average module score for cardiomyocytes
cardiomyocytes_scores <- sce@meta.data[, grep("cardiomyocytes\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_cardiomyocytes <- rowMeans(cardiomyocytes_scores)

# Calculate the average module score for epiblast
epiblast_scores <- sce@meta.data[, grep("epiblast\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_epiblast <- rowMeans(epiblast_scores)


# Flatten the list to create a vector of all module score names
all_module_score_names <- unlist(module_score_names_list)

# Create a matrix with average module scores
avg_scores <- sce@meta.data[, grep("avg_", colnames(sce@meta.data), perl = TRUE)]

# Find the dominant module for each cell
dominant_module <- colnames(avg_scores)[max.col(avg_scores)]

# Add the dominant module information to the metadata
sce@meta.data$dominant_module <- dominant_module

# Modify the dominant_module labels
dominant_module <- gsub("avg_", "", dominant_module)

# Update the metadata with the new labels
sce@meta.data$dominant_module <- dominant_module

# Choose the desired color palette
palette <- RColorBrewer::brewer.pal(7, "Dark2")  # 7 modules


plot <- DimPlot(sce, group.by = "dominant_module")
plot + scale_color_manual(values = palette) + ggtitle("Cell lineages")

# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_pluripotency", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Pluripotency")

## ******* EpCAM Add the UMAP coordinates to the metadata
sce@meta.data$`Epcam` <- sce@assays$SCT@counts["Epcam", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

epcam_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Epcam` > median(`Epcam`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Epcam expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(epcam_umap_modules)


## Dotplot heart

features <- c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog",
                       "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
                       "Fgf8", "Fgf10", "Fgf4", "Esrrb",
                       "Twist1", "Cdx2", "Hand1", "T" )

# Create the dot plot
dot_plot <- DotPlot(sce, features = features, group.by = "dominant_module") + RotatedAxis()

# Blue-to-red colour gradient
dot_plot <- dot_plot + scale_color_gradient(low = "blue", high = "red")

# Adjust plot margins / aspect ratio
dot_plot <- dot_plot + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), aspect.ratio = 0.5)

# Display the plot
print(dot_plot)


## ******* Gata6 Add the UMAP coordinates to the metadata
sce@meta.data$`Gata6` <- sce@assays$SCT@counts["Gata6", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

Gata6_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Gata6` > median(`Gata6`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Gata6 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Gata6_umap_modules)


## ******* Wnt11 Add the UMAP coordinates to the metadata
sce@meta.data$`Wnt11` <- sce@assays$SCT@counts["Wnt11", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

Wnt11_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Wnt11` > median(`Wnt11`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Wnt11 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Wnt11_umap_modules)
