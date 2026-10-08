## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 02 | Public scRNA-seq: mouse peri-implantation embryos E3.5-E6.75
##           (Mohammed et al., Cell Rep 2017; GEO: GSE100597)
## Script  : 03_pluripotency_epiblast_modules_coexpression.R
## Purpose : Final (revision) analysis: SCTransform workflow, pluripotency vs epiblast module
##           assignment, Epcam/Eomes/Foxa2/Gata6/Wnt11 expression on lineage UMAPs, dot and
##           feature plots, per-stage gene-gene Pearson correlation heatmaps and Monocle3
##           trajectory/pseudotime.
## Input   : data/Mohammed2017_GSE100597/GSE100597_count_table_QC_filtered.txt
## Output  : results/02_Mohammed2017/Mohammed2017_*.pdf, correlation_<stage>.pdf
## Notes   : choose_cells() in the Monocle3 section is interactive (opens a Shiny app).
## =============================================================================

## ---- Packages -----------------------------------------------------------------

## Cell reports cell type assigning

library(SummarizedExperiment)
library(DT)
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
library(RColorBrewer)
library(cowplot)
library(pheatmap)
library(viridis)
library(gridExtra)
library(ggridges)
library(Biobase)
library(BiocParallel)
library(purrr)
library(grid)
library(MatrixGenerics)
library(matrixStats)
library(ShortRead)
library(GenomicAlignments)
library(GenomicRanges)
library(Biostrings)
library(GenomeInfoDb)
library(XVector)
library(IRanges)
library(S4Vectors)
library(BiocGenerics)
library(Rsubread)
library(Matrix)
library(uwot)

## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/Mohammed2017_GSE100597/ (see data/README.md); all outputs are written to
## results/02_Mohammed2017/.
data_dir <- file.path(getwd(), "data", "Mohammed2017_GSE100597")
out_dir  <- file.path(getwd(), "results", "02_Mohammed2017")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Load counts, SCTransform, PCA/UMAP --------------------------------------

data <- read.table(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header = TRUE, row.names = 1, sep = "\t")
sce <- CreateSeuratObject(counts = data)

# Preprocess the data
sce <- SCTransform(sce)
sce <- FindVariableFeatures(sce)
sce <- NormalizeData(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce)

# Run UMAP
sce <- RunUMAP(sce, dims = 1:10)

dimplot <- DimPlot(sce, reduction = "umap", label = FALSE)

dimplot

# Save the combined plot
ggsave("Mohammed2017_dimplot.pdf", plot = dimplot, width = 6, height = 5, dpi = 300)


## ---- Pluripotency / epiblast module scores; dominant module per cell ------
# Add a list of known marker genes for each cell type
known_markers <- list(
  pluripotency = c("Pou5f1-17-35506037", "Sox2-3-34650005", "Nanog-6-122707489", "Lin28a-4-134003330"),
  epiblast = c("Utf1-7-139943789", "Dnmt3b-2-153649450", "L1td1-4-98726734", "Snrpn-7-59982495", "Gng3-19-8836929",
               "Gstt2-10-75831845", "Dnmt3a-12-3806007", "Dtx1-5-120680202")
)


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
palette <- RColorBrewer::brewer.pal(2, "Set1")  # 2 modules


plot <- DimPlot(sce, group.by = "dominant_module")
plot + scale_color_manual(values = palette) + ggtitle("Cell lineages")


## ---- Epcam, Foxa2, Gata6, Wnt11, Eomes on lineage UMAP (size = above median) --
## ******* EpCAM Add the UMAP coordinates to the metadata
sce@meta.data$`Epcam-17-87635979` <- sce@assays$SCT@counts["Epcam-17-87635979", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

# Create a renamed version in metadata
sce@meta.data$Epcam <- sce@meta.data$`Epcam-17-87635979`

# Re-plot using the simplified gene-name column
epcam_umap_modules <- ggplot(
  sce@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = dominant_module,
    size = Epcam > median(Epcam)
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(type = "div", palette = "Set1") +
  scale_size_manual(values = c(1.0, 2.5)) +
  labs(title = "Epcam expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(epcam_umap_modules)


# Save the combined plot
ggsave("Mohammed2017_epcam_umap_modules.pdf", plot = epcam_umap_modules, width = 6, height = 5, dpi = 300)


## ******* Foxa2 Add the UMAP coordinates to the metadata
sce@meta.data$`Foxa2-2-148042877` <- sce@assays$SCT@counts["Foxa2-2-148042877", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

# Create a renamed version in metadata
sce@meta.data$Foxa2 <- sce@meta.data$`Foxa2-2-148042877`

# Re-plot using the simplified gene-name column
Foxa2_umap_modules <- ggplot(
  sce@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = dominant_module,
    size = Foxa2 > median(Foxa2)
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(type = "div", palette = "Set1") +
  scale_size_manual(values = c(1.0, 2.5)) +
  labs(title = "Foxa2 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Foxa2_umap_modules)

# Save the combined plot
ggsave("Mohammed2017_Foxa2_umap_modules.pdf", plot = Foxa2_umap_modules, width = 6, height = 5, dpi = 300)


## ******* Gata6 Add the UMAP coordinates to the metadata
sce@meta.data$`Gata6-18-11052510` <- sce@assays$SCT@counts["Gata6-18-11052510", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

# Create a renamed version in metadata
sce@meta.data$Gata6 <- sce@meta.data$`Gata6-18-11052510`

# Re-plot using the simplified gene-name column
Gata6_umap_modules <- ggplot(
  sce@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = dominant_module,
    size = Gata6 > median(Gata6)
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(type = "div", palette = "Set1") +
  scale_size_manual(values = c(1.0, 2.5)) +
  labs(title = "Gata6 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Gata6_umap_modules)

# Save the combined plot
ggsave("Mohammed2017_Gata6_umap_modules.pdf", plot = Gata6_umap_modules, width = 6, height = 5, dpi = 300)


## ******* Wnt11 Add the UMAP coordinates to the metadata
sce@meta.data$`Wnt11-7-98835112` <- sce@assays$SCT@counts["Wnt11-7-98835112", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

# Create a renamed version in metadata
sce@meta.data$Wnt11 <- sce@meta.data$`Wnt11-7-98835112`

# Re-plot using the simplified gene-name column
Wnt11_umap_modules <- ggplot(
  sce@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = dominant_module,
    size = Wnt11 > median(Wnt11)
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(type = "div", palette = "Set1") +
  scale_size_manual(values = c(1.0, 2.5)) +
  labs(title = "Wnt11 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Wnt11_umap_modules)

# Save the combined plot
ggsave("Mohammed2017_Wnt11_umap_modules.pdf", plot = Wnt11_umap_modules, width = 6, height = 5, dpi = 300)


#### Check for Eomes gene
eomes_feat <- grep("^Eomes$|^EOMES$|Eomes|EOMES", rownames(sce@assays$SCT@counts), value = TRUE)[1]
eomes_feat


## ******* Eomes: Add expression + UMAP coords to metadata
sce@meta.data$`Eomes-9-118478212` <- sce@assays$SCT@counts["Eomes-9-118478212", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

## Create a renamed version in metadata (easier to reference in ggplot)
sce@meta.data$Eomes <- sce@meta.data$`Eomes-9-118478212`

Eomes_umap_modules <- ggplot(
  sce@meta.data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = dominant_module,
    size  = Eomes > median(Eomes)
  )
) +
  geom_point(alpha = 0.5) +
  scale_color_brewer(type = "div", palette = "Set1") +
  scale_size_manual(values = c(1.0, 2.5)) +
  labs(title = "Eomes expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

print(Eomes_umap_modules)

ggsave("Mohammed2017_Eomes_umap_modules.pdf",
       plot = Eomes_umap_modules, width = 6, height = 5, dpi = 300)


## ---- Y-chromosome gene check (Uty) ------------------------------------------
head(rownames(sce), 20)

grep("Tspx-y", rownames(sce), value = TRUE, ignore.case = TRUE)

grep("-Y-", rownames(sce), value = TRUE)


# Set SCT assay if not already set
DefaultAssay(sce) <- "SCT"

gene_names <- sapply(strsplit(rownames(sce), "-"), `[`, 1)
"Uty" %in% gene_names

# Extract gene symbols
gene_symbols <- sapply(strsplit(rownames(sce), "-"), `[`, 1)

# See if Uty or its synonym appears
"UTY" %in% gene_symbols
grep("UTY", gene_symbols, value = TRUE)


## ---- Dot plot and feature plots of EpCAM-associated genes ---------------

features <- c("sct_Epcam-17-87635979", "sct_Foxa2-2-148042877", "sct_Gata6-18-11052510", "sct_Wnt11-7-98835112",
              "sct_Nanog-6-122707489","sct_Dkk1-19-30545863", "sct_Mesp1-7-79792241", "sct_Eomes-9-118478212", "sct_Sall3-18-80966376", "sct_Mixl1-1-180693043",
              "sct_Fgf8-19-45736798", "sct_Fgf10-13-118669791", "sct_Fgf4-7-144861386", "sct_Esrrb-12-86361117",
              "sct_Twist1-12-33957671", "sct_Cdx2-5-147300805", "sct_Hand1-11-57828705", "sct_T-17-8434423")

# Create a named vector for renaming features
rename_features <- c("sct_Epcam-17-87635979" = "Epcam", "sct_Foxa2-2-148042877" = "Foxa2", "sct_Gata6-18-11052510" = "Gata6", "sct_Wnt11-7-98835112" = "Wnt11",
                     "sct_Nanog-6-122707489" = "Nanog", "sct_Dkk1-19-30545863" = "Dkk1", "sct_Mesp1-7-79792241" = "Mesp1", "sct_Eomes-9-118478212" = "Eomes", "sct_Sall3-18-80966376" = "Sall3", "sct_Mixl1-1-180693043" = "Mixl1",
                     "sct_Fgf8-19-45736798" = "Fgf8", "sct_Fgf10-13-118669791" = "Fgf10", "sct_Fgf4-7-144861386" = "Fgf4", "sct_Esrrb-12-86361117" = "Esrrb",
                     "sct_Twist1-12-33957671" = "Twist1", "sct_Cdx2-5-147300805" = "Cdx2", "sct_Hand1-11-57828705" = "Hand1", "sct_T-17-8434423" = "T")

# Create the dot plot
dot_plot <- DotPlot(sce, features = features, group.by = "dominant_module") + RotatedAxis()

# Rename the axis labels (make sure rename_features is a named vector)
dot_plot <- dot_plot + scale_x_discrete(labels = rename_features)

# Override Seurat's default color scale (this warning is normal)
dot_plot <- dot_plot + scale_color_gradient(low = "blue", high = "red")

# Fix: margin must be 4 numbers: top, right, bottom, left
dot_plot <- dot_plot + theme(
  plot.margin = unit(c(1, 1, 1, 1), "cm"),
  aspect.ratio = 0.2
)

print(dot_plot)


# Save the combined plot
ggsave("Mohammed2017_genes_dotplot.pdf", plot = dot_plot, width = 8, height = 3, dpi = 300)


## feature plot for genes
# Set SCT assay (if needed)
DefaultAssay(sce) <- "SCT"

# Loop through each feature and generate a FeaturePlot with the clean name as the title
plots <- lapply(features, function(feat) {
  FeaturePlot(sce, features = feat) +
    ggtitle(rename_features[[feat]])
})

# Show all plots in a grid
patchwork::wrap_plots(plots, ncol = 3)


# Combine all individual FeaturePlots into a grid
combined_plot <- wrap_plots(plots, ncol = 3)

# Save the combined plot
ggsave("Mohammed2017_genes_featureplot.pdf", plot = combined_plot, width = 8, height = 12, dpi = 300)


## ---- Per-stage gene-gene Pearson correlation heatmaps --------------------
####### correlation heatmaps
## Cell type
data <- read.table(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header=T, row.names=1, sep="\t")

df <-  CreateSeuratObject(counts = data, project = "Mohammed2017", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")

# Increase the time limit for subset() function
options(expressions = 10000)

# Subset cells based on quality control criteria
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 5)

df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- NormalizeData(df)

df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

## Scaling the data
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)

## Perform linear dimensional reduction
df <- RunPCA(df, features = VariableFeatures(object = df))
df <- RunUMAP(df, dims = 1:10)

# Examine and visualize PCA results a few different ways
print(df[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(df, dims = 1:2, reduction = "pca")

DimPlot(df, reduction = "umap")

cluster_list <- SplitObject(df, split.by = "orig.ident")


# Define the gene list to extract from RNA assay
genes <- c(
  "Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112",
  "Nanog-6-122707489", "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212",
  "Sall3-18-80966376", "Mixl1-1-180693043", "Fgf8-19-45736798", "Fgf10-13-118669791",
  "Fgf4-7-144861386", "Esrrb-12-86361117", "Twist1-12-33957671", "Cdx2-5-147300805",
  "Hand1-11-57828705", "T-17-8434423"
)

# Create a rename map for simplified plotting
rename_map <- setNames(
  c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog", "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
    "Fgf8", "Fgf10", "Fgf4", "Esrrb", "Twist1", "Cdx2", "Hand1", "T"),
  genes
)

# Loop through each cluster (stage) and plot correlation heatmap
for (stage in names(cluster_list)) {
  obj <- cluster_list[[stage]]

  DefaultAssay(obj) <- "RNA"

  # Filter to genes that exist in this Seurat object
  available_genes <- intersect(genes, rownames(obj))
  rename_subset <- rename_map[available_genes]

  # Skip if too few genes
  if (length(available_genes) < 2) {
    message("Skipping ", stage, ": not enough matching genes")
    next
  }

  # Get z-scored data for those genes
  expr_mat <- GetAssayData(obj, layer = "scale.data")[available_genes, ]

  # Remove genes with zero variance
  expr_mat <- expr_mat[apply(expr_mat, 1, var) > 0, , drop = FALSE]
  if (nrow(expr_mat) < 2) {
    message("Skipping ", stage, ": not enough variable genes")
    next
  }

  # Rename rows for plot labels
  rownames(expr_mat) <- rename_subset[rownames(expr_mat)]

  # Correlation matrix
  cor_mat <- cor(t(as.matrix(expr_mat)), method = "pearson")
  cor_df <- melt(cor_mat)

  # Plot
  p <- ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "gray30", mid = "white", high = "red", midpoint = 0, limits = c(-1, 1)) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.text.y = element_text(size = 10),
      panel.grid = element_blank()
    ) +
    labs(title = paste("Gene Correlation -", stage), x = NULL, y = NULL, fill = NULL)

  print(p)
  ggsave(paste0("correlation_", stage, ".pdf"), plot = p, width = 6, height = 5)
}


## ---- Monocle3 trajectory / pseudotime ------------------------------------
library(monocle3)
## Monocle analysis
# Extract counts, pData (cell data), and fData (feature data) from Seurat
counts <- sce@assays$SCT@counts
cell_data <- sce@meta.data
feature_data <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))


# Create Monocle's CellDataSet
sce_cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_data,
  gene_metadata = feature_data
)

sce_cds <- preprocess_cds(sce_cds, num_dim = 10)
sce_cds <- reduce_dimension(sce_cds, preprocess_method = 'PCA', reduction_method = 'UMAP')
sce_cds <- cluster_cells(sce_cds)
sce_cds <- learn_graph(sce_cds)
sce_cds <- order_cells(sce_cds)

names(reducedDims(sce_cds))

umap_data <- reducedDims(sce_cds)[["UMAP"]]
head(umap_data)
# Convert matrix to data frame
umap_df <- as.data.frame(umap_data)
# Merge the UMAP data frame with cell_data by row names
combined_data <- merge(umap_df, cell_data, by = "row.names", all.x = TRUE)
colnames(combined_data)[1] <- "cell_id" # Renaming the default "row.names" column


umap_data <- reducedDims(sce_cds)$UMAP
cell_data_umap <- cbind(cell_data, umap_data)


plot_pc_variance_explained(sce_cds)
plot_cells(sce_cds)

plot_cells(sce_cds, color_cells_by="orig.ident")

plot_cells(sce_cds, genes=c("Epcam-17-87635979"))

sce_cds <- reduce_dimension(sce_cds, preprocess_method = 'PCA', reduction_method="tSNE")


plot_cells(sce_cds, reduction_method="tSNE", color_cells_by="orig.ident", size=3)

p <- ggplot(umap_data, aes(x = `UMAP1`, y = `UMAP2`, color = cell_data$orig.ident)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_color_brewer(palette="Set1") +
  theme_minimal()

print(p)


plot_cells(sce_cds, color_cells_by="orig.ident", label_groups_by_cluster=TRUE)


## Marker gene expression
marker_test_res <- top_markers(sce_cds, group_cells_by="orig.ident",
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

unique(colData(sce_cds)$partition)
table(top_specific_markers$cell_group)

## top genes in each stages
plot_genes_by_group(sce_cds,
                    top_specific_marker_ids,
                    group_cells_by="orig.ident",
                    ordering_type="maximal_on_diag",
                    max.size=3)

## top specific marker
top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

plot_genes_by_group(sce_cds,
                    top_specific_marker_ids,
                    group_cells_by="orig.ident",
                    ordering_type="cluster_row_col",
                    max.size=3)


## subset
sce_cds_subset <- choose_cells(sce_cds)


plot_cells(sce_cds)
plot_cells(sce_cds, color_cells_by = "pseudotime")

sce_cds <- order_cells(sce_cds, reduction_method = "UMAP")

sce_cds <- learn_graph(sce_cds, use_partition = TRUE)


plot_cells(sce_cds, reduction_method = 'UMAP', color_cells_by = "pseudotime")

plot_cells(sce_cds, color_cells_by = "orig.ident")

plot_cells(sce_cds, show_trajectory_graph = TRUE, color_cells_by = "pseudotime")

plot_cells(sce_cds, genes = c('Epcam-17-87635979'))


cell_data <- colData(sce_cds)
ggplot(cell_data, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point() +
  scale_color_viridis_c()


colnames(colData(sce_cds))
unique(colData(sce_cds)$orig.ident)


plot_cells(sce_cds, color_cells_by = "orig.ident")


## Publication-style plot of embryonic stages
# Plot cells
p <- plot_cells(sce_cds, color_cells_by = "orig.ident")

# Adjust theme elements
p <- p + theme_minimal() +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  labs(title = "Embryonic stages")

print(p)


plot_cells(sce_cds, show_trajectory_graph = TRUE, color_cells_by = "pseudotime")


## ---- Refined module scoring and publication dot plot ---------------------
# ------------------- Define known marker genes -------------------
known_markers <- list(
  pluripotency = c("Pou5f1-17-35506037", "Sox2-3-34650005", "Nanog-6-122707489", "Lin28a-4-134003330"),
  epiblast = c("Utf1-7-139943789", "Dnmt3b-2-153649450", "L1td1-4-98726734", "Snrpn-7-59982495", "Gng3-19-8836929",
               "Gstt2-10-75831845", "Dnmt3a-12-3806007", "Dtx1-5-120680202")
)

# ------------------- Calculate module scores -------------------
for (cell_type in names(known_markers)) {
  sce <- AddModuleScore(sce, features = known_markers[[cell_type]], name = paste0(cell_type, "_"))
}

# ------------------- Compute average scores -------------------
sce$avg_pluripotency <- rowMeans(sce@meta.data[, grep("pluripotency_\\d+", colnames(sce@meta.data))])
sce$avg_epiblast     <- rowMeans(sce@meta.data[, grep("epiblast_\\d+", colnames(sce@meta.data))])

# ------------------- Assign dominant module -------------------
score_matrix <- sce@meta.data[, c("avg_pluripotency", "avg_epiblast")]
sce$dominant_module <- gsub("avg_", "", colnames(score_matrix)[max.col(score_matrix)])

# ------------------- Visualize cell type identity -------------------
palette <- RColorBrewer::brewer.pal(2, "Set1")
DimPlot(sce, group.by = "dominant_module") + scale_color_manual(values = palette) + ggtitle("Cell lineages")

# ------------------- Gene features for DotPlot -------------------
features <- c(
  "Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112",
  "Nanog-6-122707489", "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212",
  "Sall3-18-80966376", "Mixl1-1-180693043", "Fgf8-19-45736798", "Fgf10-13-118669791",
  "Fgf4-7-144861386", "Esrrb-12-86361117", "Twist1-12-33957671", "Cdx2-5-147300805",
  "Hand1-11-57828705", "T-17-8434423"
)

rename_features <- setNames(
  c("Epcam", "Foxa2", "Gata6", "Wnt11", "Nanog", "Dkk1", "Mesp1", "Eomes", "Sall3", "Mixl1",
    "Fgf8", "Fgf10", "Fgf4", "Esrrb", "Twist1", "Cdx2", "Hand1", "T"),
  features
)

# ------------------- Remove duplicate gene columns from meta.data -------------------
genes_to_remove <- intersect(features, colnames(sce@meta.data))
sce@meta.data[, genes_to_remove] <- NULL


dot_data <- DotPlot(sce, features = features, group.by = "dominant_module")$data
head(dot_data)

subset(dot_data, features.plot == "Epcam-17-87635979")

subset(dot_data, features.plot == "Epcam")


# Use raw dot plot data
dot_data <- DotPlot(sce, features = features, group.by = "dominant_module")$data
dot_data$features.plot <- rename_features[as.character(dot_data$features.plot)]

# Custom dot plot with styling
pub_dotplot <- ggplot(dot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp)) +
  scale_size(range = c(2, 8), name = "Percent Expressed") +
  scale_color_gradient(low = "blue", high = "red", name = "Avg. Expression") +
  labs(x = "Features", y = "Identity") +
  theme_minimal(base_size = 14) +
  theme(
    # Axis text
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),

    # Axis titles
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),

    # Remove grid
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # Add axis lines
    axis.line = element_line(color = "black", linewidth = 0.5),

    # Legends and spacing
    legend.title = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 13),
    legend.position = "right",
    plot.margin = unit(c(1, 1, 1, 1), "cm"),

    # Optional aspect ratio
    aspect.ratio = 0.25
  )

# Display the plot
print(pub_dotplot)


# Save the combined plot
ggsave("Mohammed2017_genes_dotplot_publication.pdf", plot = pub_dotplot, width = 11, height = 5, dpi = 300)


# Check what assay and slot are used
DefaultAssay(sce)
# Check what layer DotPlot is pulling from (uses data slot by default)
head(sce@assays$SCT@data["Epcam-17-87635979", ])

DefaultAssay(sce) <- "SCT"
VlnPlot(sce, features = "Epcam-17-87635979", group.by = "dominant_module", assay = "SCT")

summary(sce@meta.data$avg_pluripotency)
summary(sce@meta.data$avg_epiblast)
table(sce@meta.data$dominant_module)

FeatureScatter(sce, feature1 = "avg_pluripotency", feature2 = "avg_epiblast") +
  ggtitle("Module score comparison")


## feature plot for genes
# Set SCT assay (if needed)
DefaultAssay(sce) <- "SCT"

# Loop through each feature and generate a FeaturePlot with the clean name as the title
plots <- lapply(features, function(feat) {
  FeaturePlot(sce, features = feat) +
    ggtitle(rename_features[[feat]])
})

# Show all plots in a grid
patchwork::wrap_plots(plots, ncol = 3)


# Combine all individual FeaturePlots into a grid
combined_plot <- wrap_plots(plots, ncol = 3)


# Save the combined plot
ggsave("Mohammed2017_genes_featureplot_refined.pdf", plot = combined_plot, width = 8, height = 12, dpi = 300)
