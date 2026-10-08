## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 02 | Public scRNA-seq: mouse peri-implantation embryos E3.5-E6.75
##           (Mohammed et al., Cell Rep 2017; GEO: GSE100597)
## Script  : 01_Seurat_QC_marker_gene_expression.R
## Purpose : Seurat (LogNormalize) QC, normalisation, PCA/UMAP and expression of the EpCAM-
##           associated genes (Epcam, Eomes, Foxa2, Gata6, Wnt11, Nanog, ...) across E3.5-E6.75;
##           export of normalised expression; clustering and cluster markers.
## Input   : data/Mohammed2017_GSE100597/GSE100597_count_table_QC_filtered.txt
## Output  : results/02_Mohammed2017/ (normalized_expression.*, expression_data.csv, Seurat .rds)
## Notes   : Gene identifiers in GSE100597 are 'Symbol-chr-start' (e.g. Epcam-17-87635979).
## =============================================================================

## ---- Packages -----------------------------------------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/Mohammed2017_GSE100597/ (see data/README.md); all outputs are written to
## results/02_Mohammed2017/.
data_dir <- file.path(getwd(), "data", "Mohammed2017_GSE100597")
out_dir  <- file.path(getwd(), "results", "02_Mohammed2017")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Load counts, QC filtering, LogNormalize, PCA/UMAP ----------------------

data <- read.table(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header=T, row.names=1, sep="\t")

df <-  CreateSeuratObject(counts = data, project = "Mohammed2017", min.cells = 3, min.features = 200)

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
write.table(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expression.txt", sep = "\t", quote = FALSE)

# Save normalized expression data as a CSV file
write.csv(x = as.matrix(GetAssayData(object = df, slot = "data")), file = "normalized_expression.csv", row.names = TRUE)

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

df <- RunUMAP(df, dims = 1:10)


## ---- Stage marker genes (identities = embryonic stage) -------------------
markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

## ---- Expression of EpCAM-associated genes (PCA feature, ridge, violin, dot plots) --
## Training codes
## Select gene of interest
gene_of_interest <- "Epcam-17-87635979"

## Visualize gene expression on PCA plot
FeaturePlot(df, features = gene_of_interest, reduction = "pca")

## Nanog
genes_of_interest <- c("Nanog-6-122707489")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Gata6
genes_of_interest <- c("Gata6-18-11052510")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Foxa2
genes_of_interest <- c("Foxa2-2-148042877")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Wnt11
genes_of_interest <- c("Wnt11-7-98835112")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Dkk1
genes_of_interest <- c("Dkk1-19-30545863")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Mesp1
genes_of_interest <- c("Mesp1-7-79792241")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Eomes
genes_of_interest <- c("Eomes-9-118478212")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Sall3
genes_of_interest <- c("Sall3-18-80966376")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Mixl1
genes_of_interest <- c("Mixl1-1-180693043")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Foxc1
genes_of_interest <- c("Foxc1-13-31806646")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Stat3
genes_of_interest <- c("Stat3-11-100885098")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## fgf8
genes_of_interest <- c("Fgf8-19-45736798")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## fgf10
genes_of_interest <- c("Fgf10-13-118669791")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Klf5
genes_of_interest <- c("Klf5-14-99298691")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Col5a2
genes_of_interest <- c("Col5a2-1-45374321")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Col5a1
genes_of_interest <- c("Col5a1-2-27886425")
FeaturePlot(df, features = genes_of_interest, reduction = "pca")

## Overall marker genes
genes_of_interest <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112", "Nanog-6-122707489",
                       "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
                       "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Fgf4-7-144861386", "Esrrb-12-86361117",
                      "Twist1-12-33957671", "Cdx2-5-147300805", "Hand1-11-57828705", "T-17-8434423" )

FeaturePlot(df, features = genes_of_interest, reduction = "pca")


## Ridge plot selected genes
features <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112")
RidgePlot(df, features = features, ncol = 2)
VlnPlot(df, features = features)

## Expression data for the gene of interest
genes_of_interest <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112", "Nanog-6-122707489",
                       "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
                       "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Fgf4-7-144861386", "Esrrb-12-86361117",
                       "Twist1-12-33957671", "Cdx2-5-147300805", "Hand1-11-57828705", "T-17-8434423")

expression_data <- as.data.frame(GetAssayData(object = df, slot = "data")[genes_of_interest, ])

write.csv(expression_data, file = "expression_data.csv")


features <- c("Epcam-17-87635979", "Nanog-6-122707489", "Gata6-18-11052510", "Foxa2-2-148042877")
df
# Ridge plots - from ggridges. Visualize single cell expression distributions in each cluster
RidgePlot(df, features = features, ncol = 2)

features <- c("Wnt11-7-98835112", "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212")
RidgePlot(df, features = features, ncol = 2)

features <- c("Sall3-18-80966376", "Mixl1-1-180693043",  "Foxc1-13-31806646", "Stat3-11-100885098")
RidgePlot(df, features = features, ncol = 2)

features <- c("Fgf8-19-45736798", "Fgf10-13-118669791", "Klf5-14-99298691", "Col5a2-1-45374321", "Col5a1-2-27886425")
RidgePlot(df, features = features, ncol = 2)

## Violin plots
features <- c("Epcam-17-87635979", "Nanog-6-122707489", "Gata6-18-11052510", "Foxa2-2-148042877", "Wnt11-7-98835112",
              "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
              "Foxc1-13-31806646", "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Klf5-14-99298691",
              "Col5a2-1-45374321", "Col5a1-2-27886425")
# Violin plot - Visualize single cell expression distributions in each cluster
VlnPlot(df, features = features)

# Dot plots - the size of the dot corresponds to the percentage of cells expressing the feature
# in each cluster. The color represents the average expression level
DotPlot(df, features = features) + RotatedAxis()

# Single cell heatmap of feature expression
DoHeatmap(subset(df, downsample = 100), features = features, size = 3)

## ---- Dimensionality, clustering and cluster markers ------------------------
## Dimension plots
DimHeatmap(df, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(df, dims = 1:15, cells = 500, balanced = TRUE)

## Determine the ‘dimensionality’ of the dataset
df <- JackStraw(df, num.replicate = 100)
df <- ScoreJackStraw(df, dims = 1:20)

JackStrawPlot(df, dims = 1:15)

ElbowPlot(df)

## Cluster the cells
df <- FindNeighbors(df, dims = 1:10)
df <- FindClusters(df, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(df), 5)

df <- RunUMAP(df, dims = 1:10)

# individual clusters
DimPlot(df, reduction = "umap")

## Save RDS
saveRDS(df, file = "Mohammed2017_seurat_LogNormalize.rds")

## Finding differentially expressed features (cluster biomarkers)
# find all markers of cluster 2
cluster2.markers <- FindMarkers(df, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(df, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
df.markers <- FindAllMarkers(df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
df.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

cluster0.markers <- FindMarkers(df, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

VlnPlot(df, features = c("Epcam-17-87635979", "Gata6-18-11052510"))
