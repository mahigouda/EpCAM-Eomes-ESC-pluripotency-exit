# =============================================================================
# Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
#           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
# Module  : 04 | Public scRNA-seq: human embryos, Carnegie stages CS12-CS16
#           (Xu et al., Nat Cell Biol 2023; GEO: GSE157329)
# Script  : 01_human_embryo_marker_coexpression.py
# Purpose : Scanpy analysis of human embryo scRNA-seq (CS12-CS16): QC, normalisation, UMAP by
#           stage / developmental system; expression of EPCAM and EpCAM-associated genes
#           (EOMES, FOXA2, GATA6, WNT11, NANOG, ...), cell counts per stage and system, and
#           single-cell co-expression of gene combinations (e.g. EPCAM+EOMES+...) per stage.
# Input   : data/Xu2023_GSE157329/{GSE157329_raw_counts.mtx, GSE157329_gene_annotate.txt,
#           #           barcode_new_1.txt, GSE157329_cell_annotate_1.txt}
# Output  : results/04_human_embryo_GSE157329/ (figures/*.pdf, *.pdf, *_merged_data_stages.csv)
# Notes   : Converted from the original Jupyter notebook; run cell-by-cell or as a script.
# =============================================================================

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import os
from scipy.io import mmread
import anndata
import anndata as ad
import scipy.io
import igraph
from scipy.stats import pearsonr
import louvain
import leidenalg
import seaborn as sns
import matplotlib.colors as mcolors

# ---- Paths --------------------------------------------------------------------
# Run from the repository root. Inputs in data/Xu2023_GSE157329/ (see data/README.md);
# outputs (figures, tables) are written to results/04_human_embryo_GSE157329/.
from pathlib import Path
DATA_DIR = Path("data/Xu2023_GSE157329").resolve()
OUT_DIR = Path("results/04_human_embryo_GSE157329").resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(OUT_DIR)
sc.settings.figdir = "figures"  # scanpy `save=` figures go to results/04_human_embryo_GSE157329/figures/

adata = sc.read_mtx(DATA_DIR / "GSE157329_raw_counts.mtx")
adata = adata.T

adata_1 = adata.copy()

# Read the gene information file
genes_df = pd.read_csv(DATA_DIR / "GSE157329_gene_annotate.txt", header=None, sep='\t', index_col=0,
                       names=['gene_ids', 'gene_symbols'])
genes_df.set_index('gene_symbols', inplace=True)
genes_df = genes_df.iloc[1:]
genes_df.index.name = 'gene_ids'

# Assign gene_ids using the index column of genes_df
adata_1.var['gene_ids'] = genes_df.index

barcodes = pd.read_csv(DATA_DIR / "barcode_new_1.txt", header=None, sep='\t', names=['barcode'])
# Remove the first row in the DataFrame
barcodes = barcodes.iloc[1:]

adata_1.obs.index = barcodes['barcode']

metadata = pd.read_csv(DATA_DIR / "GSE157329_cell_annotate_1.txt", sep='\t')
metadata = metadata.set_index("barcode")

common_rows = adata_1.obs.index.intersection(metadata.index)

columns_to_add = ['cell_id', 'cluster_id', 'developmental system', 'annotation',
          'final_annotation', 'embryo', 'sample', 'stage', 'dissection_part',
          'total_UMIs']

adata_1.obs.loc[common_rows, columns_to_add] = metadata.loc[common_rows, columns_to_add]

adata_2 = adata_1.copy()

sc.pp.filter_cells(adata_2, min_genes = 200)
sc.pp.filter_cells(adata_2, max_genes=5000)
sc.pp.filter_genes(adata_2, min_cells=3)
sc.pp.normalize_total(adata_2, target_sum=1e4)
sc.pp.log1p(adata_2)
sc.pp.highly_variable_genes(adata_2, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata_2, max_value=10)
sc.tl.pca(adata_2, svd_solver='arpack')
sc.pp.neighbors(adata_2)
sc.tl.louvain(adata_2)
sc.tl.umap(adata_2)

# Generate the UMAP plot
plt.rcParams.update({'font.size': 18})
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.linewidth'] = 0

sc.pl.umap(adata_2, color=['stage'], show=True, size=5, save="_Nature_human_embryo_scRNA_stages.pdf")
sc.pl.umap(adata_2, color=['dissection_part'], show=True, size=5, save="_Nature_human_embryo_scRNA_dissection_part.pdf")
sc.pl.umap(adata_2, color=['sample'], show=True, size=5, save="_Nature_human_embryo_scRNA_sample.pdf")
sc.pl.umap(adata_2, color=['developmental system'], show=True, size=5, save="_Nature_human_embryo_scRNA_developmental system.pdf")
sc.pl.umap(adata_2, color=['embryo'], show=True, size=5, save="_Nature_human_embryo_scRNA_embryo.pdf")

# Access the gene_ids index from adata_2.var
gene_ids = adata_2.var.index
print(gene_ids[:100])

# Set the gene symbols as the index of the .var DataFrame
adata_2.var.index = adata_2.var['gene_ids']

# Access the gene_ids index from adata_2.var
gene_ids = adata_2.var.index
print(gene_ids[:100])

genes_of_interest = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG",
                     "DKK1", "MESP1", "EOMES", "SALL3", "MIXL1",
                     "STAT3", "FGF8", "FGF10", "FGF4", "ESRRB",
                     "TWIST1", "CDX2", "HAND1", "TBXT"]

# Keep only genes present in the dataset
genes_of_interest_present = [gene for gene in genes_of_interest if gene in adata_2.var.index]

# Define the custom color map
colors = ["#cccccc", "#ff6666", "#ff0000"]
cmap = mcolors.LinearSegmentedColormap.from_list("", colors)

# Plot the UMAP with the custom color map
sc.pl.umap(adata_2, color=genes_of_interest_present, cmap=cmap, size=5, save="_Nature_human_embryo_scRNA_genes_of_interest.pdf")

print(adata_2.obs)

adata_2_genes_of_interest = adata_2[:, adata_2.var_names.isin(genes_of_interest_present)]

cell_counts = adata_2_genes_of_interest.obs.groupby(['stage']).size()
print(cell_counts)

# Subset adata to only include cells that express 'Epcam'
adata_2_Epcam = adata_2[adata_2[:, 'EPCAM'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Epcam = adata_2_Epcam.obs.groupby(['stage']).size()
print(cell_counts_Epcam)


# Subset adata to only include cells that express 'Foxa2'
adata_2_Foxa2 = adata_2[adata_2[:, 'FOXA2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Foxa2 = adata_2_Foxa2.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Gat6'
adata_2_Gata6 = adata_2[adata_2[:, 'GATA6'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Gata6 = adata_2_Gata6.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'WNT11'
adata_2_Wnt11 = adata_2[adata_2[:, 'WNT11'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Wnt11 = adata_2_Wnt11.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Nanog'
adata_2_Nanog = adata_2[adata_2[:, 'NANOG'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Nanog = adata_2_Nanog.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'DKK1'
adata_2_Dkk1 = adata_2[adata_2[:, 'DKK1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Dkk1 = adata_2_Dkk1.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Mesp1'
adata_2_Mesp1 = adata_2[adata_2[:, 'MESP1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Mesp1 = adata_2_Mesp1.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Eomes'
adata_2_Eomes = adata_2[adata_2[:, 'EOMES'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Eomes = adata_2_Eomes.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Sall3'
adata_2_Sall3 = adata_2[adata_2[:, 'SALL3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Sall3 = adata_2_Sall3.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Mixl1'
adata_2_Mixl1 = adata_2[adata_2[:, 'MIXL1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Mixl1 = adata_2_Mixl1.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Stat3'
adata_2_Stat3 = adata_2[adata_2[:, 'STAT3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Stat3 = adata_2_Stat3.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Fgf8'
adata_2_Fgf8 = adata_2[adata_2[:, 'FGF8'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Fgf8 = adata_2_Fgf8.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Fgf10'
adata_2_Fgf10 = adata_2[adata_2[:, 'FGF10'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Fgf10 = adata_2_Fgf10.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Fgf4'
adata_2_Fgf4 = adata_2[adata_2[:, 'FGF4'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Fgf4 = adata_2_Fgf4.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Esrrb'
adata_2_Esrrb = adata_2[adata_2[:, 'ESRRB'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Esrrb = adata_2_Esrrb.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Twist1'
adata_2_Twist1 = adata_2[adata_2[:, 'TWIST1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Twist1 = adata_2_Twist1.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Cdx2'
adata_2_Cdx2 = adata_2[adata_2[:, 'CDX2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Cdx2 = adata_2_Cdx2.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'Hand1'
adata_2_Hand1 = adata_2[adata_2[:, 'HAND1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Hand1 = adata_2_Hand1.obs.groupby(['stage']).size()

# Subset adata to only include cells that express 'TBXT'
adata_2_Tbxt = adata_2[adata_2[:, 'TBXT'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Tbxt = adata_2_Tbxt.obs.groupby(['stage']).size()

# Combine all cell counts into one DataFrame
stages = pd.concat([cell_counts_Epcam, cell_counts_Foxa2, cell_counts_Gata6, cell_counts_Wnt11, cell_counts_Nanog, cell_counts_Dkk1,
                   cell_counts_Mesp1, cell_counts_Eomes, cell_counts_Sall3, cell_counts_Mixl1,
                   cell_counts_Stat3, cell_counts_Fgf8, cell_counts_Fgf10, cell_counts_Fgf4, cell_counts_Esrrb,
                   cell_counts_Twist1, cell_counts_Cdx2, cell_counts_Hand1, cell_counts_Tbxt], axis=1)
stages.columns = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG",
                     "DKK1", "MESP1", "EOMES", "SALL3", "MIXL1",
                     "STAT3", "FGF8", "FGF10", "FGF4", "ESRRB",
                     "TWIST1", "CDX2", "HAND1", "TBXT"]

# Melt the DataFrame to long format for easier plotting
stages_melted = stages.reset_index().melt(id_vars='stage', var_name='gene', value_name='cell_count')

# Relabel the stages
stages_melted['stage'] = stages_melted['stage'].replace({
    'CS13-14': 'CS13_14',
    'CS15-16': 'CS15_16'
})

# Define the gene groups
genes_group1 = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG", "DKK1",]
genes_group2 = ["MESP1", "EOMES", "SALL3", "MIXL1", "STAT3", "FGF8",]
genes_group3 = ["FGF10", "FGF4", "ESRRB", "TWIST1", "CDX2", "HAND1", "TBXT"]

# Define stage order manually
stage_order = ['CS12', 'CS13_14', 'CS15_16']


# Convert 'stage' to an ordered categorical type
stages_melted['stage'] = pd.Categorical(stages_melted['stage'], categories=stage_order, ordered=True)

# Convert 'stage' to numerical values for smoothing
stages_melted['stage_num'] = stages_melted['stage'].str.extract(r'(\d+.\d+)').astype(float)

# Sort data for rolling mean
stages_melted.sort_values(by=['gene', 'stage_num'], inplace=True)

# Compute rolling mean with a window of 3
stages_melted['cell_count_smooth'] = stages_melted.groupby('gene')['cell_count'].transform(lambda x: x.rolling(3, 1).mean())


## Line plot verticle

fig, axs = plt.subplots(3, figsize=(6, 12))

for ax, genes_group in zip(axs.flatten(), [genes_group1, genes_group2, genes_group3]):
    # Filter the data for the current gene group
    stages_melted_group = stages_melted[stages_melted['gene'].isin(genes_group)]

    # Create the line plot for the current gene group with smoothed cell counts
    sns.lineplot(x='stage', y='cell_count_smooth', hue='gene', data=stages_melted_group, ax=ax)

    # Add a title and labels
    ax.set_ylabel('Cell Counts', fontsize=16)
    ax.set_xlabel('Stage', fontsize=16)

    # Rotate the x-axis labels for better readability and increase label size
    ax.tick_params(axis='x', rotation=90, labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # Increase the size of the legend, add a title, and move it to the right of the plot
    ax.legend(title='Gene', title_fontsize='14', fontsize='14', bbox_to_anchor=(1.05, 1), loc='upper left')

# Remove the top and right spines for a cleaner look
sns.despine()

# Ensure the layout is well-arranged
plt.tight_layout()

# Save the plot as a PDF
plt.savefig('Human_embryo_Cell_count_plot_Stages.pdf', format='pdf')

# Show the plot
plt.show()


## Line plot Horizontal
# Set the number of rows and columns
num_rows = 1
num_cols = 3

# Adjust the figsize based on the desired arrangement
fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, 6))

# Flatten the axes array if necessary
axs = axs.flatten()

for i, genes_group in enumerate([genes_group1, genes_group2, genes_group3]):
    # Filter the data for the current gene group
    stages_melted_group = stages_melted[stages_melted['gene'].isin(genes_group)]

    # Create the line plot for the current gene group with smoothed cell counts
    sns.lineplot(x='stage', y='cell_count_smooth', hue='gene', data=stages_melted_group, ax=axs[i])

    # Add a title and labels
    axs[i].set_ylabel('Cell Counts', fontsize=16)
    axs[i].set_xlabel('Stage', fontsize=16)

    # Rotate the x-axis labels for better readability and increase label size
    axs[i].tick_params(axis='x', rotation=90, labelsize=14)
    axs[i].tick_params(axis='y', labelsize=14)

    # Increase the size of the legend, add a title, and move it to the right of the plot
    axs[i].legend(title='Gene', title_fontsize='14', fontsize='14', bbox_to_anchor=(1.05, 1), loc='upper left')

# Remove the top and right spines for a cleaner look
sns.despine()

# Ensure the layout is well-arranged
plt.tight_layout()

# Save the plot as a PDF
plt.savefig('Human_embryo_Cell_count_plot_Stages_Horizontal.pdf', format='pdf')

# Show the plot
plt.show()


## Cell counts developmental system
cell_counts_developmental_system = adata_2_genes_of_interest.obs.groupby(['developmental system']).size()
print(cell_counts_developmental_system)

# Subset adata to only include cells that express 'Epcam'
adata_celltype_Epcam = adata_2[adata_2[:, 'EPCAM'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Epcam = adata_celltype_Epcam.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Foxa2'
adata_celltype_Foxa2 = adata_2[adata_2[:, 'FOXA2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Foxa2 = adata_celltype_Foxa2.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Gata6'
adata_celltype_Gata6 = adata_2[adata_2[:, 'GATA6'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Gata6 = adata_celltype_Gata6.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Wnt11'
adata_celltype_Wnt11 = adata_2[adata_2[:, 'WNT11'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Wnt11 = adata_celltype_Wnt11.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Nanog'
adata_celltype_Nanog = adata_2[adata_2[:, 'NANOG'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Nanog = adata_celltype_Nanog.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Dkk1'
adata_celltype_Dkk1 = adata_2[adata_2[:, 'DKK1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Dkk1 = adata_celltype_Dkk1.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Mesp1'
adata_celltype_Mesp1 = adata_2[adata_2[:, 'MESP1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Mesp1 = adata_celltype_Mesp1.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Eomes'
adata_celltype_Eomes = adata_2[adata_2[:, 'EOMES'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Eomes = adata_celltype_Eomes.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Sall3'
adata_celltype_Sall3 = adata_2[adata_2[:, 'SALL3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Sall3 = adata_celltype_Sall3.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Mixl1'
adata_celltype_Mixl1 = adata_2[adata_2[:, 'MIXL1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Mixl1 = adata_celltype_Mixl1.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Stat3'
adata_celltype_Stat3 = adata_2[adata_2[:, 'STAT3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Stat3 = adata_celltype_Stat3.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Fgf8'
adata_celltype_Fgf8 = adata_2[adata_2[:, 'FGF8'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Fgf8 = adata_celltype_Fgf8.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Fgf10'
adata_celltype_Fgf10 = adata_2[adata_2[:, 'FGF10'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Fgf10 = adata_celltype_Fgf10.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Fgf4'
adata_celltype_Fgf4 = adata_2[adata_2[:, 'FGF4'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Fgf4 = adata_celltype_Fgf4.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Esrrb'
adata_celltype_Esrrb = adata_2[adata_2[:, 'ESRRB'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Esrrb = adata_celltype_Esrrb.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Twist1'
adata_celltype_Twist1 = adata_2[adata_2[:, 'TWIST1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Twist1 = adata_celltype_Twist1.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Cdx2'
adata_celltype_Cdx2 = adata_2[adata_2[:, 'CDX2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Cdx2 = adata_celltype_Cdx2.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'Hand1'
adata_celltype_Hand1 = adata_2[adata_2[:, 'HAND1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Hand1 = adata_celltype_Hand1.obs.groupby(['developmental system']).size()

# Subset adata to only include cells that express 'T'
adata_celltype_Tbxt = adata_2[adata_2[:, 'TBXT'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_celltype_Tbxt = adata_celltype_Tbxt.obs.groupby(['developmental system']).size()


# Combine all cell counts into one DataFrame
celltypes = pd.concat([cell_counts_celltype_Epcam, cell_counts_celltype_Foxa2, cell_counts_celltype_Gata6,
                       cell_counts_celltype_Wnt11, cell_counts_celltype_Nanog, cell_counts_celltype_Dkk1,
                   cell_counts_celltype_Mesp1, cell_counts_celltype_Eomes, cell_counts_celltype_Sall3,
                       cell_counts_celltype_Mixl1, cell_counts_celltype_Stat3, cell_counts_celltype_Fgf8,
                       cell_counts_celltype_Fgf10, cell_counts_celltype_Fgf4, cell_counts_celltype_Esrrb,
                   cell_counts_celltype_Twist1, cell_counts_celltype_Cdx2, cell_counts_celltype_Hand1,
                       cell_counts_celltype_Tbxt], axis=1)
celltypes.columns = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG",
                     "DKK1", "MESP1", "EOMES", "SALL3", "MIXL1",
                     "STAT3", "FGF8", "FGF10", "FGF4", "ESRRB",
                     "TWIST1", "CDX2", "HAND1", "TBXT"]

# Melt the DataFrame to long format for easier plotting
celltypes_melted = celltypes.reset_index().melt(id_vars='developmental system', var_name='gene', value_name='cell_count')

celltypes_melted.head


## Dot plot
# Define size factor for your dots. You might have to adjust this value to get your desired output.
size_factor = 0.3

# Create a scatter plot
plt.figure(figsize=(42, 30))
scatter_plot = sns.scatterplot(x='gene', y='developmental system', size='cell_count',
                               sizes=(0.3, celltypes_melted['cell_count'].max() * size_factor),
                               data=celltypes_melted, color = 'red')

# Increase the size of the axis labels
plt.tick_params(axis='both', which='major', labelsize=32)

# Remove x-axis and y-axis labels
plt.xlabel('')
plt.ylabel('')

# Access the legend from the scatter plot, set label size and move it to the right outside of the plot
legend = scatter_plot.legend(loc='center right', bbox_to_anchor=(1.07, 0.5), ncol=1,
                             labelspacing=1.5, handletextpad=1.5)  # Adjust these values as needed
for text in legend.texts:
    text.set_fontsize(32)

# Save the plot as a PDF
plt.savefig('Human_embryo_dotplot_celltypes_cell_counts.pdf', format='pdf')

# Show the plot
plt.show()


## Coexpression
signature_genes = ["EPCAM", "FOXA2", "GATA6", "WNT11"]
# Increase font size for axis labels and colorbar
plt.rcParams.update({'font.size': 18})
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.linewidth'] = 0

# Define the custom color map
colors = ["#cccccc", "#ff6666", "#ff0000"]
cmap = mcolors.LinearSegmentedColormap.from_list("", colors)

# Plot the UMAP with the custom color map
sc.pl.umap(adata_2, color=signature_genes, cmap=cmap, size=5, save="_Human_embryo_Epcam_Foxa2_Gata6_Wnt11_.pdf")


# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks = [adata_2[:, gene].X > 0 for gene in signature_genes]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask = np.all(masks, axis=0)

# Convert the boolean mask to integer
adata_2.obs['coexpressed'] = coexpressed_mask.astype(int)

# Subset adata to only include cells that coexpress all the signature genes
adata_2_coexpressed = adata_2[adata_2.obs['coexpressed'] == 1]

# Count the number of cells in each stage
cell_counts_coexpressed = adata_2_coexpressed.obs.groupby(['stage']).size()
print(cell_counts_coexpressed)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Epcam_Foxa2_Gata6_Wnt11.pdf")


Epcam_Foxa2 = ["EPCAM", "FOXA2"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2 = [adata_2[:, gene].X > 0 for gene in Epcam_Foxa2]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2 = np.all(masks_Epcam_Foxa2, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Foxa2'] = coexpressed_mask_Epcam_Foxa2.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Foxa2', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Epcam_Foxa2.pdf")


Epcam_Gata6 = ["EPCAM", "GATA6"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Gata6 = [adata_2[:, gene].X > 0 for gene in Epcam_Gata6]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Gata6 = np.all(masks_Epcam_Gata6, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Gata6'] = coexpressed_mask_Epcam_Gata6.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Gata6', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Epcam_Gata6.pdf")


Epcam_Wnt11 = ["EPCAM", "WNT11"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Wnt11 = [adata_2[:, gene].X > 0 for gene in Epcam_Wnt11]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Wnt11 = np.all(masks_Epcam_Wnt11, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Wnt11'] = coexpressed_mask_Epcam_Wnt11.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Wnt11', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Epcam_Wnt11.pdf")


Foxa2_Gata6 = ["FOXA2", "GATA6"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Foxa2_Gata6 = [adata_2[:, gene].X > 0 for gene in Foxa2_Gata6]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Foxa2_Gata6 = np.all(masks_Foxa2_Gata6, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Foxa2_Gata6'] = coexpressed_mask_Foxa2_Gata6.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Foxa2_Gata6', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Foxa2_Gata6.pdf")


Foxa2_Wnt11 = ["FOXA2", "WNT11"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Foxa2_Wnt11 = [adata_2[:, gene].X > 0 for gene in Foxa2_Wnt11]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Foxa2_Wnt11 = np.all(masks_Foxa2_Wnt11, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Foxa2_Wnt11'] = coexpressed_mask_Foxa2_Wnt11.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Foxa2_Wnt11', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Foxa2_Wnt11.pdf")


Gata6_Wnt11 = ["GATA6", "WNT11"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Gata6_Wnt11 = [adata_2[:, gene].X > 0 for gene in Gata6_Wnt11]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Gata6_Wnt11 = np.all(masks_Gata6_Wnt11, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Gata6_Wnt11'] = coexpressed_mask_Gata6_Wnt11.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Gata6_Wnt11', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Gata6_Wnt11.pdf")


Epcam_Foxa2_Gata6 = ["EPCAM", "FOXA2", "GATA6"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2_Gata6 = [adata_2[:, gene].X > 0 for gene in Epcam_Foxa2_Gata6]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2_Gata6 = np.all(masks_Epcam_Foxa2_Gata6, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Foxa2_Gata6'] = coexpressed_mask_Epcam_Foxa2_Gata6.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Foxa2_Gata6', cmap=cmap, size=5, save="_Human_embryo_coexpressed_Epcam_Foxa2_Gata6.pdf")


Epcam_Nanog_Eomes_Fgf8 = ["EPCAM", "NANOG", "EOMES", "FGF8"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Nanog_Eomes_Fgf8 = [adata_2[:, gene].X > 0 for gene in Epcam_Nanog_Eomes_Fgf8]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Nanog_Eomes_Fgf8 = np.all(masks_Epcam_Nanog_Eomes_Fgf8, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Nanog_Eomes_Fgf8'] = coexpressed_mask_Epcam_Nanog_Eomes_Fgf8.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
# UMAP with fixed colour scale (0/1)
sc.pl.umap(adata_2, color='coexpressed_Epcam_Nanog_Eomes_Fgf8', cmap=cmap, size=5,
           vmin=0, vmax=1,
           save="_Human_embryo_coexpressed_Epcam_Nanog_Eomes_Fgf8.pdf")

Epcam_Eomes_Cdx2_Hand1 = ["EPCAM", "EOMES", "CDX2", "HAND1"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Eomes_Cdx2_Hand1 = [adata_2[:, gene].X > 0 for gene in Epcam_Eomes_Cdx2_Hand1]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Eomes_Cdx2_Hand1 = np.all(masks_Epcam_Eomes_Cdx2_Hand1, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Eomes_Cdx2_Hand1'] = coexpressed_mask_Epcam_Eomes_Cdx2_Hand1.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Eomes_Cdx2_Hand1', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Eomes_Cdx2_Hand1.pdf")


Epcam_Foxa2_Gata6_Stat3_Cdx2 = ["EPCAM", "FOXA2", "GATA6", "STAT3", "CDX2"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2_Gata6_Stat3_Cdx2 = [adata_2[:, gene].X > 0 for gene in Epcam_Foxa2_Gata6_Stat3_Cdx2]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2_Gata6_Stat3_Cdx2 = np.all(masks_Epcam_Foxa2_Gata6_Stat3_Cdx2, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Foxa2_Gata6_Stat3_Cdx2'] = coexpressed_mask_Epcam_Foxa2_Gata6_Stat3_Cdx2.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Foxa2_Gata6_Stat3_Cdx2', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Foxa2_Gata6_Stat3_Cdx2.pdf")


Epcam_Foxa2_Stat3_Cdx2 = ["EPCAM", "FOXA2", "STAT3", "CDX2"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2_Stat3_Cdx2 = [adata_2[:, gene].X > 0 for gene in Epcam_Foxa2_Stat3_Cdx2]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2_Stat3_Cdx2 = np.all(masks_Epcam_Foxa2_Stat3_Cdx2, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Foxa2_Stat3_Cdx2'] = coexpressed_mask_Epcam_Foxa2_Stat3_Cdx2.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Foxa2_Stat3_Cdx2', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Foxa2_Stat3_Cdx2.pdf")


Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T = ["EPCAM", "GATA6", "MESP1", "EOMES", "MIXL1", "FGF8", "TBXT"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T = [adata_2[:, gene].X > 0 for gene in Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T = np.all(masks_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T'] = coexpressed_mask_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_T.pdf")


Epcam_Nanog_Eomes_Mixl1_Fgf8_T = ["EPCAM", "NANOG", "EOMES", "MIXL1", "FGF8", "TBXT"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Nanog_Eomes_Mixl1_Fgf8_T = [adata_2[:, gene].X > 0 for gene in Epcam_Nanog_Eomes_Mixl1_Fgf8_T]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Nanog_Eomes_Mixl1_Fgf8_T = np.all(masks_Epcam_Nanog_Eomes_Mixl1_Fgf8_T, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Nanog_Eomes_Mixl1_Fgf8_T'] = coexpressed_mask_Epcam_Nanog_Eomes_Mixl1_Fgf8_T.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Nanog_Eomes_Mixl1_Fgf8_T', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Nanog_Eomes_Mixl1_Fgf8_T.pdf")


Epcam_Dkk1_Stat3_Ffg8 = ["EPCAM", "DKK1", "STAT3", "FGF8"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Dkk1_Stat3_Ffg8 = [adata_2[:, gene].X > 0 for gene in Epcam_Dkk1_Stat3_Ffg8]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Dkk1_Stat3_Ffg8 = np.all(masks_Epcam_Dkk1_Stat3_Ffg8, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Dkk1_Stat3_Ffg8'] = coexpressed_mask_Epcam_Dkk1_Stat3_Ffg8.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Dkk1_Stat3_Ffg8', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Dkk1_Stat3_Ffg8.pdf")


Epcam_Dkk1 = ["EPCAM", "DKK1"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Dkk1 = [adata_2[:, gene].X > 0 for gene in Epcam_Dkk1]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Dkk1 = np.all(masks_Epcam_Dkk1, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Dkk1'] = coexpressed_mask_Epcam_Dkk1.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Dkk1', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Dkk1.pdf")


Epcam_Stat3 = ["EPCAM", "STAT3"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Stat3 = [adata_2[:, gene].X > 0 for gene in Epcam_Stat3]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Stat3 = np.all(masks_Epcam_Stat3, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Stat3'] = coexpressed_mask_Epcam_Stat3.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Stat3', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Stat3.pdf")


Epcam_Fgf8 = ["EPCAM", "FGF8"]
# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Fgf8 = [adata_2[:, gene].X > 0 for gene in Epcam_Fgf8]
# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Fgf8 = np.all(masks_Epcam_Fgf8, axis=0)
# Convert the boolean mask to integer
adata_2.obs['coexpressed_Epcam_Fgf8'] = coexpressed_mask_Epcam_Fgf8.astype(int)
# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_2, color='coexpressed_Epcam_Fgf8', cmap=cmap, size=5,
           save="_Human_embryo_coexpressed_Epcam_Fgf8.pdf")


# Subset adata to only include cells that coexpress "Epcam", "Foxa2", "Gata6", and "Wnt11"
adata_celltype_Epcam_Foxa2_Gata6_Wnt11 = adata_2[
    np.logical_and.reduce(
        (
            adata_2[:, "EPCAM"].X > 0,
            adata_2[:, "FOXA2"].X > 0,
           adata_2[:, "GATA6"].X > 0,
            adata_2[:, "WNT11"].X > 0,
        )
    )
]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2_Gata6_Wnt11 = adata_celltype_Epcam_Foxa2_Gata6_Wnt11.obs.groupby(["stage"]).size()


# Subset adata to only include cells that coexpress "Epcam" and "Foxa2"
adata_celltype_Epcam_Foxa2 = adata_2[np.logical_and(adata_2[:, "EPCAM"].X > 0,
                                                           adata_2[:, "FOXA2"].X > 0)]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2 = adata_celltype_Epcam_Foxa2.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Gata6"
adata_celltype_Epcam_Gata6 = adata_2[np.logical_and(adata_2[:, "EPCAM"].X > 0, adata_2[:, "GATA6"].X > 0)]
# Count the number of cells in each stage
cell_counts_Epcam_Gata6 = adata_celltype_Epcam_Gata6.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Foxa2" and "Gata6"
adata_celltype_Foxa2_Gata6 = adata_2[np.logical_and(adata_2[:, "FOXA2"].X > 0, adata_2[:, "GATA6"].X > 0)]
# Count the number of cells in each stage
cell_counts_Foxa2_Gata6 = adata_celltype_Foxa2_Gata6.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Foxa2", and "Gata6"
adata_celltype_Epcam_Foxa2_Gata6 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                  adata_2[:, "FOXA2"].X > 0,
                                                                  adata_2[:, "GATA6"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2_Gata6 = adata_celltype_Epcam_Foxa2_Gata6.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Nanog", "Eomes" and "Fgf8"
adata_celltype_Epcam_Nanog_Eomes_Fgf8 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "NANOG"].X > 0,
                                                                         adata_2[:, "EOMES"].X > 0,
                                                                         adata_2[:, "FGF8"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Nanog_Eomes_Fgf8 = adata_celltype_Epcam_Nanog_Eomes_Fgf8.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Eomes", "Cdx2", and "Hand1"
adata_celltype_Epcam_Eomes_Cdx2_Hand1 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "EOMES"].X > 0,
                                                                         adata_2[:, "CDX2"].X > 0,
                                                                         adata_2[:, "HAND1"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Eomes_Cdx2_Hand1 = adata_celltype_Epcam_Eomes_Cdx2_Hand1.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Foxa2", "Stat3", and "Cdx2"
adata_celltype_Epcam_Foxa2_Stat3_Cdx2 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "FOXA2"].X > 0,
                                                                         adata_2[:, "STAT3"].X > 0,
                                                                         adata_2[:, "CDX2"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2_Stat3_Cdx2 = adata_celltype_Epcam_Foxa2_Stat3_Cdx2.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Gata6", "Mesp1", "Eomes", "Mixl1", "Fgf8", and "Tbxt"
adata_celltype_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_Tbxt = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "GATA6"].X > 0,
                                                                         adata_2[:, "MESP1"].X > 0,
                                                                         adata_2[:, "EOMES"].X > 0,
                                                                         adata_2[:, "MIXL1"].X > 0,
                                                                         adata_2[:, "FGF8"].X > 0,
                                                                         adata_2[:, "TBXT"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_Tbxt = adata_celltype_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_Tbxt.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Nanog", "Eomes", "Mixl1", "Fgf8", and "Tbxt"
adata_celltype_Epcam_Nanog_Eomes_Mixl1_Fgf8_Tbxt = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "NANOG"].X > 0,
                                                                         adata_2[:, "EOMES"].X > 0,
                                                                         adata_2[:, "MIXL1"].X > 0,
                                                                         adata_2[:, "FGF8"].X > 0,
                                                                         adata_2[:, "TBXT"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Nanog_Eomes_Mixl1_Fgf8_Tbxt = adata_celltype_Epcam_Nanog_Eomes_Mixl1_Fgf8_Tbxt.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Dkk1", "Stat3", and "Fgf8"
adata_celltype_Epcam_Dkk1_Stat3_Ffg8 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0,
                                                                         adata_2[:, "DKK1"].X > 0,
                                                                         adata_2[:, "STAT3"].X > 0,
                                                                         adata_2[:, "FGF8"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Dkk1_Stat3_Ffg8 = adata_celltype_Epcam_Dkk1_Stat3_Ffg8.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Dkk1"
adata_celltype_Epcam_Dkk1 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0, adata_2[:, "DKK1"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Dkk1 = adata_celltype_Epcam_Dkk1.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Stat3"
adata_celltype_Epcam_Stat3 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0, adata_2[:, "STAT3"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Stat3 = adata_celltype_Epcam_Stat3.obs.groupby(['stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Fgf8"
adata_celltype_Epcam_Fgf8 = adata_2[np.logical_and.reduce((adata_2[:, "EPCAM"].X > 0, adata_2[:, "FGF8"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Fgf8 = adata_celltype_Epcam_Fgf8.obs.groupby(['stage']).size()


## Coexpression data merging
merged_data = pd.concat([cell_counts_Epcam_Foxa2,
                         cell_counts_Epcam_Gata6, cell_counts_Foxa2_Gata6,
                         cell_counts_Epcam_Foxa2_Gata6, cell_counts_Epcam_Nanog_Eomes_Fgf8,
                         cell_counts_Epcam_Eomes_Cdx2_Hand1, cell_counts_Epcam_Foxa2_Stat3_Cdx2,
                         cell_counts_Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_Tbxt,
                         cell_counts_Epcam_Nanog_Eomes_Mixl1_Fgf8_Tbxt, cell_counts_Epcam_Dkk1_Stat3_Ffg8,
                         cell_counts_Epcam_Dkk1, cell_counts_Epcam_Stat3, cell_counts_Epcam_Fgf8], axis=1)

merged_data.columns = ['Epcam_Foxa2', 'Epcam_Gata6', 'Foxa2_Gata6', 'Epcam_Foxa2_Gata6', 'Epcam_Nanog_Eomes_Fgf8',
                       'Epcam_Eomes_Cdx2_Hand1', 'Epcam_Foxa2_Stat3_Cdx2', 'Epcam_Gata6_Mesp1_Eomes_Mixl1_Fgf8_Tbxt',
                       'Epcam_Nanog_Eomes_Mixl1_Fgf8_Tbxt', 'Epcam_Dkk1_Stat3_Ffg8', 'Epcam_Dkk1', 'Epcam_Stat3', 'Epcam_Fgf8']

merged_data.fillna(0, inplace=True)
merged_data.to_csv('Human_embryo_gene_coexpression_merged_data_stages.csv', index=False)


# Relabel the stages in the index
merged_data = merged_data.rename(index={
    'CS13-14': 'CS13_14',
    'CS15-16': 'CS15_16'
})

stage_order = ['CS12', 'CS13_14', 'CS15_16']

## Lineplot for coexpressed single cells
# Ensure that 'stage' is of type 'category' and is ordered
merged_data.index = pd.CategoricalIndex(merged_data.index, categories=stage_order, ordered=True)

# Sort the dataframe by 'stage'
merged_data = merged_data.sort_index()

# Create a figure and axes
fig, ax = plt.subplots(figsize=(14, 6))

# Define the palette to ensure each gene has a distinct color
palette = sns.color_palette("husl", n_colors=len(merged_data.columns) - 1)

# Iterate over each gene column (excluding the 'stage' column) and plot the line
for column, color in zip(merged_data.columns[1:], palette):
    sns.lineplot(data=merged_data, x=merged_data.index, y=column, marker='o', ax=ax, color=color, label=column)

# Add a title and labels
ax.set_title('Coexpressed single cells', fontsize=20)
ax.set_ylabel('Cell Counts', fontsize=16)
ax.set_xlabel('Stage', fontsize=16)

# Rotate the x-axis labels for better readability and increase label size
ax.tick_params(axis='x', rotation=90, labelsize=14)
ax.tick_params(axis='y', labelsize=14)

# Increase the size of the legend and add a title
ax.legend(title='Gene', title_fontsize='14', fontsize='14', bbox_to_anchor=(1.05, 1), loc='upper left')

# Remove the top and right spines for a cleaner look
sns.despine()

# Ensure the layout is well-arranged
plt.tight_layout()

# Save the plot as a PDF
plt.savefig('Human_embryo_Lineplot_Coexpression_singlecells_count.pdf', format='pdf')

# Show the plot
plt.show()


## Heatmap for coexpressed single cells
# Ensure that 'stage' is of type 'category' and is ordered
merged_data.index = pd.CategoricalIndex(merged_data.index, categories=stage_order, ordered=True)

# Sort the dataframe by 'stage'
merged_data = merged_data.sort_index()

# Create a figure and axes
fig, ax = plt.subplots(figsize=(20, 12))

# Define the custom color map
colors = ["#cccccc", "#ff6666", "#ff0000"]
cmap = mcolors.LinearSegmentedColormap.from_list("", colors)

# Generate the heatmap
heatmap = sns.heatmap(merged_data, cmap=cmap, annot=True, fmt='.2f', ax=ax)

# Set the title and labels
ax.set_title('Coexpressed single cells', fontsize=20)
ax.set_xlabel('Genes', fontsize=16)
ax.set_ylabel('Stage', fontsize=16)

# Rotate the x-axis labels for better readability
plt.xticks(rotation=45, ha='right')

# Show the colorbar
cbar = heatmap.collections[0].colorbar
cbar.set_label('Cell Counts', fontsize=14)

# Adjust the layout
plt.tight_layout()

# Save the plot as a PDF
plt.savefig('Human_embryo_coexpression_new_heatmap.pdf', format='pdf')

# Show the plot
plt.show()
