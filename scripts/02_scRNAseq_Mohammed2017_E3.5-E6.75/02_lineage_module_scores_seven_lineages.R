## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 02 | Public scRNA-seq: mouse peri-implantation embryos E3.5-E6.75
##           (Mohammed et al., Cell Rep 2017; GEO: GSE100597)
## Script  : 02_lineage_module_scores_seven_lineages.R
## Purpose : Assign each cell to a dominant lineage module (pluripotency, epiblast, mesoendoderm,
##           mesoderm, endoderm, ectoderm, cardiomyocytes) with Seurat::AddModuleScore and
##           visualise Epcam and EpCAM-associated genes across lineages (UMAP, dot plots).
## Input   : data/Mohammed2017_GSE100597/GSE100597_count_table_QC_filtered.txt
## Output  : results/02_Mohammed2017/dot_plot_data.csv (+ interactive plots)
## =============================================================================

## ---- Packages -----------------------------------------------------------------

## Cell reports cell type assigning

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(patchwork)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(scater)
library(reshape2)
library(ggplot2)
library(monocle3)
library(RColorBrewer)
library(cowplot)
library(pheatmap)


## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/Mohammed2017_GSE100597/ (see data/README.md); all outputs are written to
## results/02_Mohammed2017/.
data_dir <- file.path(getwd(), "data", "Mohammed2017_GSE100597")
out_dir  <- file.path(getwd(), "results", "02_Mohammed2017")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Load counts, SCTransform, PCA/UMAP --------------------------------------

# Load the data
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

# Visualize the results
DimPlot(sce, reduction = "umap")

genes_of_interest <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112", "Nanog-6-122707489",
                       "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
                       "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Fgf4-7-144861386", "Esrrb-12-86361117",
                       "Twist1-12-33957671", "Cdx2-5-147300805", "Hand1-11-57828705", "T-17-8434423" )

FeaturePlot(sce, features = genes_of_interest, reduction = "umap")


## ---- Module scores for seven lineages; dominant module per cell ---------

# Add a list of known marker genes for each cell type
known_markers <- list(
  pluripotency = c("Pou5f1-17-35506037", "Sox2-3-34650005", "Nanog-6-122707489", "Lin28a-4-134003330"),
  mesoendoderm = c("Mixl1-1-180693043", "Gsc-12-104471209", "Foxa2-2-148042877"),
  mesoderm = c("Tbx6-7-126781483", "Mesp1-7-79792241", "Hand1-11-57828705", "Myf5-10-107482908"),
  endoderm = c("Sox17-1-4490928", "Foxa2-2-148042877", "Sox7-14-63943674", "Gata4-14-63198922"),
  ectoderm = c("Pax6-2-105668896", "Neurod1-2-79452521", "Nes-3-87971078"),
  cardiomyocytes = c("Myl4-11-104550663", "Tnnt2-1-135836386", "Myl7-11-5896637", "Tnni1-1-135783065", "Nexn-3-152236986",
                     "Nkx2-5-17-26838664", "Acta2-19-34241090", "Mef2c-13-83504034", "Actc1-2-114047282", "Gyg-3-20122084"),
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

# Calculate the average module score for mesoendoderm
mesoendoderm_scores <- sce@meta.data[, grep("mesoendoderm\\d+", colnames(sce@meta.data), perl = TRUE)]
sce@meta.data$avg_mesoendoderm <- rowMeans(mesoendoderm_scores)

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

## ******* EpCAM Add the UMAP coordinates to the metadata
sce@meta.data$`Epcam-17-87635979` <- sce@assays$SCT@counts["Epcam-17-87635979", ]

sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

epcam_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Epcam-17-87635979` > median(`Epcam-17-87635979`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Epcam expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(epcam_umap_modules)


## ---- Dot plot of EpCAM-associated genes per lineage ----------------------
## Dot plot Epcam
# Set the features you want to visualize

features <- c("sct_Epcam-17-87635979", "sct_Foxa2-2-148042877", "sct_Gata6-18-11052510", "sct_Wnt11-7-98835112",
              "sct_Nanog-6-122707489","sct_Dkk1-19-30545863", "sct_Mesp1-7-79792241", "sct_Eomes-9-118478212", "sct_Sall3-18-80966376", "sct_Mixl1-1-180693043",
              "sct_Stat3-11-100885098", "sct_Fgf8-19-45736798", "sct_Fgf10-13-118669791", "sct_Fgf4-7-144861386", "sct_Esrrb-12-86361117",
              "sct_Twist1-12-33957671", "sct_Cdx2-5-147300805", "sct_Hand1-11-57828705", "sct_T-17-8434423")

# Create a named vector for renaming features
rename_features <- c("sct_Epcam-17-87635979" = "Epcam", "sct_Foxa2-2-148042877" = "Foxa2", "sct_Gata6-18-11052510" = "Gata6", "sct_Wnt11-7-98835112" = "Wnt11",
                     "sct_Nanog-6-122707489" = "Nanog", "sct_Dkk1-19-30545863" = "Dkk1", "sct_Mesp1-7-79792241" = "Mesp1", "sct_Eomes-9-118478212" = "Eomes", "sct_Sall3-18-80966376" = "Sall3", "sct_Mixl1-1-180693043" = "Mixl1",
                     "sct_Stat3-11-100885098" = "Stat3", "sct_Fgf8-19-45736798" = "Fgf8", "sct_Fgf10-13-118669791" = "Fgf10", "sct_Fgf4-7-144861386" = "Fgf4", "sct_Esrrb-12-86361117" = "Esrrb",
                     "sct_Twist1-12-33957671" = "Twist1", "sct_Cdx2-5-147300805" = "Cdx2", "sct_Hand1-11-57828705" = "Hand1", "sct_T-17-8434423" = "T")


# Create the dot plot
dot_plot <- DotPlot(sce, features = features, group.by = "dominant_module") + RotatedAxis()

# Rename the axis labels
dot_plot <- dot_plot + scale_x_discrete(labels = rename_features)

# Blue-to-red colour gradient
dot_plot <- dot_plot + scale_color_gradient(low = "blue", high = "red")

# Adjust plot margins / aspect ratio
dot_plot <- dot_plot + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), aspect.ratio = 0.3)

# Display the plot
print(dot_plot)


# Filter the cells based on the dominant_module column
sce_filtered <- sce[, which(sce$dominant_module %in% c("pluripotency", "epiblast"))]

# Create the dot plot
dot_plot <- DotPlot(sce_filtered, features = features, group.by = "dominant_module") + RotatedAxis()

# Rename the axis labels
dot_plot <- dot_plot + scale_x_discrete(labels = rename_features)

# Blue-to-red colour gradient
dot_plot <- dot_plot + scale_color_gradient(low = "blue", high = "red")

# Adjust plot margins / aspect ratio
dot_plot <- dot_plot + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"), aspect.ratio = 0.2)

# Display the plot
print(dot_plot)

# Save the dot_plot data as a CSV file
write.csv(dot_plot$data, file = "dot_plot_data.csv")


## ---- Module-score feature plots and lineage UMAPs ---------------------

# Visualize the module scores on the UMAP plot
FeaturePlot(sce, features = all_module_score_names, pt.size = 0.5, cols = c("lightgrey", "blue", "red"))


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_pluripotency", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Pluripotency")

# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_mesoendoderm", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Mesoendoderm")


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_mesoderm", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Mesoderm")


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_endoderm", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Endoderm")


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_ectoderm", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Ectoderm")


# Create a matrix with average module scores
avg_scores <- sce@meta.data[, grep("avg_", colnames(sce@meta.data), perl = TRUE)]

# Find the dominant module for each cell
dominant_module <- colnames(avg_scores)[max.col(avg_scores)]

# Add the dominant module information to the metadata
sce@meta.data$dominant_module <- dominant_module

# Visualize all modules in one UMAP plot
DimPlot(sce, group.by = "dominant_module")


# Modify the dominant_module labels
dominant_module <- gsub("avg_", "", dominant_module)

# Update the metadata with the new labels
sce@meta.data$dominant_module <- dominant_module


# Visualize all modules in one UMAP plot with updated labels
DimPlot(sce, group.by = "dominant_module")


# Choose the desired color palette
palette <- RColorBrewer::brewer.pal(5, "Set1")  # 5 modules


plot <- DimPlot(sce, group.by = "dominant_module")
plot + scale_color_manual(values = palette) + ggtitle("Cell lineages")


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_cardiomyocytes", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Cardiomyocytes")


# Visualize the average module score on the UMAP plot
plot <- FeaturePlot(sce, features = "avg_epiblast", pt.size = 0.5, cols = c("lightgrey", "grey", "magenta"))

# Add a title to the plot
plot + ggtitle("Epiblast")


# Create the UMAP plot of cell lineages
umap_lineages <- DimPlot(sce, group.by = "dominant_module") +
  scale_color_manual(values = palette) +
  ggtitle("Cell lineages")

# Create the UMAP plot of Epcam-17-87635979 gene expression
epcam_umap <- FeaturePlot(sce, features = "sct_Epcam-17-87635979", reduction = "umap", pt.size = 0.1) +
  ggtitle("Epcam-17-87635979 expression")

# Combine the UMAP plot of cell lineages and Epcam-17-87635979 expression
combined_umap <- cowplot::plot_grid(umap_lineages, epcam_umap, ncol = 2)

# Display the combined UMAP plot
print(combined_umap)

umap_lineages <- DimPlot(sce, group.by = "dominant_module", label = TRUE) +
  scale_color_manual(values = palette) +
  ggtitle("Cell lineages and Epcam-17-87635979 expression")

epcam_expression <- FeaturePlot(sce, features = "sct_Epcam-17-87635979", reduction = "umap", pt.size = 0.1)

# Combine the plots
merged_plot <- umap_lineages + epcam_expression + plot_layout(guides = "collect")

# Display the merged plot
print(merged_plot)


# Create a new SingleCellExperiment object with Epcam-17-87635979 expression as a metadata column
sce_with_epcam <- sce
sce_with_epcam@meta.data$Epcam <- sce_with_epcam@assays$RNA@data["Epcam-17-87635979",]


## ---- Gene expression (size = above-median) on lineage UMAP -------------


colnames(sce@meta.data)
head(sce@meta.data)

## ******* EpCAM Add the UMAP coordinates to the metadata
sce@meta.data$UMAP_1 <- sce@reductions$umap@cell.embeddings[, 1]
sce@meta.data$UMAP_2 <- sce@reductions$umap@cell.embeddings[, 2]

epcam_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Epcam-17-87635979` > median(`Epcam-17-87635979`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Epcam expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(epcam_umap_modules)


## *******Foxa2 Add the UMAP coordinates to the metadata
Foxa2_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Foxa2-2-148042877` > median(`Foxa2-2-148042877`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Foxa2 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Foxa2_umap_modules)

## *******Gata6 Add the UMAP coordinates to the metadata
Gata6_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Gata6-18-11052510` > median(`Gata6-18-11052510`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Gata6 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Gata6_umap_modules)

## *******Wnt11 Add the UMAP coordinates to the metadata
Wnt11_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Wnt11-7-98835112` > median(`Wnt11-7-98835112`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Wnt11 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Wnt11_umap_modules)

## *******Nanog Add the UMAP coordinates to the metadata
Nanog_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Nanog-6-122707489` > median(`Nanog-6-122707489`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Nanog expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Nanog_umap_modules)

## *******Dkk1 Add the UMAP coordinates to the metadata
Dkk1_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Dkk1-19-30545863` > median(`Dkk1-19-30545863`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Dkk1 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Dkk1_umap_modules)


## ******* Mesp1 Add the UMAP coordinates to the metadata
Mesp1_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Mesp1-7-79792241` > median(`Mesp1-7-79792241`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Mesp1 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Mesp1_umap_modules)

## ******* Eomes Add the UMAP coordinates to the metadata
Eomes_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Eomes-9-118478212` > median(`Eomes-9-118478212`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Eomes expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Eomes_umap_modules)

## ******* Sall3 Add the UMAP coordinates to the metadata
Sall3_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Sall3-18-80966376` > median(`Sall3-18-80966376`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Sall3 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Sall3_umap_modules)

## ******* Mixl1 Add the UMAP coordinates to the metadata
Mixl1_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Mixl1-1-180693043` > median(`Mixl1-1-180693043`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Mixl1 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Mixl1_umap_modules)

## ******* Stat3 Add the UMAP coordinates to the metadata
Stat3_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Stat3-11-100885098` > median(`Stat3-11-100885098`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Stat3 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Stat3_umap_modules)

## ******* Fgf8 Add the UMAP coordinates to the metadata
Fgf8_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Fgf8-19-45736798` > median(`Fgf8-19-45736798`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Fgf8 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Fgf8_umap_modules)

## ******* Fgf10 Add the UMAP coordinates to the metadata
Fgf10_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Fgf10-13-118669791` > median(`Fgf10-13-118669791`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Fgf10 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Fgf10_umap_modules)

## ******* Fgf4 Add the UMAP coordinates to the metadata
Fgf4_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Fgf4-7-144861386` > median(`Fgf4-7-144861386`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Fgf4 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Fgf4_umap_modules)

## ******* Esrrb Add the UMAP coordinates to the metadata
Esrrb_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Esrrb-12-86361117` > median(`Esrrb-12-86361117`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Esrrb expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Esrrb_umap_modules)

## ******* Twist1 Add the UMAP coordinates to the metadata
Twist1_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Twist1-12-33957671` > median(`Twist1-12-33957671`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Twist1 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Twist1_umap_modules)

## ******* Cdx2 Add the UMAP coordinates to the metadata
Cdx2_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Cdx2-5-147300805` > median(`Cdx2-5-147300805`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Cdx2 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Cdx2_umap_modules)

## ******* Hand1 Add the UMAP coordinates to the metadata
Hand1_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `Hand1-11-57828705` > median(`Hand1-11-57828705`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Hand1 expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(Hand1_umap_modules)

## ******* T Add the UMAP coordinates to the metadata
T_umap_modules <- ggplot(sce@meta.data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = `T-17-8434423` > median(`T-17-8434423`))) +
  geom_point(alpha = 0.8) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "T expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(T_umap_modules)


## ******** Marker genes expression in Cell lineages
merged_plot <- plot_grid(
  epcam_umap_modules, Gata6_umap_modules, Wnt11_umap_modules, Nanog_umap_modules, Dkk1_umap_modules,
  Mesp1_umap_modules, Eomes_umap_modules, Sall3_umap_modules, Mixl1_umap_modules, Stat3_umap_modules,
  Fgf8_umap_modules, Fgf10_umap_modules, Fgf4_umap_modules, Esrrb_umap_modules, Twist1_umap_modules,
  Cdx2_umap_modules, Hand1_umap_modules, T_umap_modules, ncol = 5, nrow = 4
)

# Display the merged plot
print(merged_plot)


genes <- c("Epcam-17-87635979", "Foxa2-2-148042877", "Gata6-18-11052510", "Wnt11-7-98835112", "Nanog-6-122707489",
           "Dkk1-19-30545863", "Mesp1-7-79792241", "Eomes-9-118478212", "Sall3-18-80966376", "Mixl1-1-180693043",
           "Stat3-11-100885098", "Fgf8-19-45736798", "Fgf10-13-118669791", "Fgf4-7-144861386", "Esrrb-12-86361117",
           "Twist1-12-33957671", "Cdx2-5-147300805", "Hand1-11-57828705", "T-17-8434423" )


library(gtable)

# Prepare the data for plotting
plot_data <- sce@meta.data %>%
  dplyr::select(UMAP_1, UMAP_2, dominant_module, `Epcam-17-87635979`, `Foxa2-2-148042877`, `Gata6-18-11052510`, `Wnt11-7-98835112`, `Nanog-6-122707489`, `Dkk1-19-30545863`, `Mesp1-7-79792241`, `Eomes-9-118478212`, `Sall3-18-80966376`, `Mixl1-1-180693043`, `Stat3-11-100885098`, `Fgf8-19-45736798`, `Fgf10-13-118669791`, `Fgf4-7-144861386`, `Esrrb-12-86361117`, `Twist1-12-33957671`, `Cdx2-5-147300805`, `Hand1-11-57828705`, `T-17-8434423`) %>%
  tidyr::gather(gene, expression, -UMAP_1, -UMAP_2, -dominant_module)

# Create the custom UMAP plot
umap_plot <- ggplot(plot_data, aes(x = UMAP_1, y = UMAP_2, color = dominant_module, size = expression > median(expression))) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ gene) +
  scale_size_manual(values = c(1, 2.5)) +
  labs(title = "Gene expression in Cell lineages") +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(umap_plot)


## ---- LogNormalize workflow and marker availability check ---------------

data <- read.table(file.path(data_dir, "GSE100597_count_table_QC_filtered.txt"), header=T, row.names=1, sep="\t")

df <-  CreateSeuratObject(counts = data, project = "Cell_reports", min.cells = 3, min.features = 200)

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")

# Increase the time limit for subset() function
options(expressions = 10000)

# Subset cells based on quality control criteria
df <- subset(df, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5)

df <- NormalizeData(df, normalization.method = "LogNormalize", scale.factor = 10000)
df <- NormalizeData(df)

df <- FindVariableFeatures(df, selection.method = "vst", nfeatures = 2000)

## Scaling the data
all.genes <- rownames(df)
df <- ScaleData(df, features = all.genes)

## Perform linear dimensional reduction
df <- RunPCA(df, features = VariableFeatures(object = df))

# Examine and visualize PCA results a few different ways
print(df[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(df, dims = 1:2, reduction = "pca")

DimPlot(df, reduction = "pca")

# subset the Seurat object to keep only cells in the E3.5 cluster
df_E3.5 <- subset(df, idents = "E3.5")


# Identify marker genes for each lineage
pluripotency_markers <- c("Pou5f1", "Sox2", "Nanog", "Lin28a")
mesoendoderm_markers <- c("Mixl1", "Gsc", "Foxa2", "Sox17")
mesoderm_markers <- c("Tbx6", "Mesp1", "Hand1", "Myf5")
endoderm_markers <- c("Sox17", "Foxa2", "Sox7", "Gata4")
ectoderm_markers <- c("Sox1", "Pax6", "Neurod1", "Nes")

genes <- rownames(df)
pluripotency_markers[!pluripotency_markers %in% genes]
mesoendoderm_markers[!mesoendoderm_markers %in% genes]
mesoderm_markers[!mesoderm_markers %in% genes]
endoderm_markers[!endoderm_markers %in% genes]
ectoderm_markers[!ectoderm_markers %in% genes]


## ---- Per-stage objects (QC: nFeature < 7500) --------------------------
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

# Examine and visualize PCA results a few different ways
print(df[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(df, dims = 1:2, reduction = "pca")

DimPlot(df, reduction = "pca")

cluster_list <- SplitObject(df, split.by = "orig.ident")
