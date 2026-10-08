# =============================================================================
# Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
#           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
# Module  : 05 | Public scRNA-seq: zebrafish embryogenesis (Drop-seq, URD lineage atlas;
#           Farrell et al., Science 2018; GEO: GSE106587)
# Script  : 01_zebrafish_marker_coexpression.py
# Purpose : Scanpy analysis of the zebrafish URD Drop-seq atlas: lineage annotation from URD
#           lineage flags, UMAP by lineage/stage, expression of epcam and EpCAM-associated
#           orthologues, cell counts per stage/lineage and single-cell co-expression per stage.
# Input   : data/Farrell2018_zebrafish_URD/{URD_Dropseq_Expression_Log2TPM.txt, meta.txt}
# Output  : results/05_zebrafish/ (figures/*.pdf, Zebrafish_*.pdf, Zebrafish_coexpression_merged_data_stages.csv)
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
import seaborn as sns
from matplotlib.cm import get_cmap
import gzip
import tarfile
import louvain
import leidenalg
import matplotlib.colors as mcolors
from scipy.sparse import issparse

# ---- Paths --------------------------------------------------------------------
# Run from the repository root. Inputs in data/Farrell2018_zebrafish_URD/ (see data/README.md);
# outputs (figures, tables) are written to results/05_zebrafish/.
from pathlib import Path
DATA_DIR = Path("data/Farrell2018_zebrafish_URD").resolve()
OUT_DIR = Path("results/05_zebrafish").resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(OUT_DIR)
sc.settings.figdir = "figures"  # scanpy `save=` figures go to results/05_zebrafish/figures/


## Load expression matrix (log2 TPM) and cell metadata
data_new = pd.read_csv(DATA_DIR / "URD_Dropseq_Expression_Log2TPM.txt", header=None, sep='\t', dtype=str)
# Set the first row as column names
data_new.columns = data_new.iloc[0]

# Drop the first row as it's now the column names
data_new = data_new.iloc[1:]

# Set an appropriate column as the index (e.g., using the 'NAME' column)
data_new = data_new.set_index('GENE')

data_new = data_new.T

metadata = pd.read_csv(DATA_DIR / "meta.txt", header=None, sep='\t', dtype=str)
print(metadata.head())
# Set the first row as column names
metadata.columns = metadata.iloc[0]
# Drop the first row as it's now the column names
metadata = metadata.iloc[1:]
# Set an appropriate column as the index (e.g., using the 'NAME' column)
metadata = metadata.set_index('NAME')


adata = sc.AnnData(data_new)
# Update the adata.obs DataFrame with the annotations
adata.obs = adata.obs.merge(metadata, left_index=True, right_index=True)


adata_1 = adata.copy()
# Assuming adata_1.X is a sparse matrix
if issparse(adata_1.X):
    adata_1.X = adata_1.X.astype(np.float64)

# Convert string values to float, while handling non-numeric entries
adata_1.X = np.array(adata_1.X, dtype=np.float64)

# Define the list of lineage columns


lineage_columns = [
   'Lineage_Spinal_Cord', 'Lineage_Diencephalon', 'Lineage_Optic_Cup', 'Lineage_Midbrain_Neural_Crest', 'Lineage_Hindbrain_R3',
    'Lineage_Hindbrain_R4_5_6', 'Lineage_Telencephalon', 'Lineage_Epidermis', 'Lineage_Neural_Plate_Border', 'Lineage_Placode_Adeno._Lens_Trigeminal',
    'Lineage_Placode_Epibranchial_Otic', 'Lineage_Placode_Olfactory', 'Lineage_Tailbud', 'Lineage_Adaxial_Cells', 'Lineage_Somites', 'Lineage_Hematopoeitic_ICM',
    'Lineage_Hematopoeitic_RBI_Pronephros', 'Lineage_Endoderm_Pharyngeal', 'Lineage_Endoderm_Pancreatic_Intestinal', 'Lineage_Heart_Primordium',
    'Lineage_Cephalic_Mesoderm', 'Lineage_Prechordal_Plate', 'Lineage_Notochord', 'Lineage_Primordial_Germ_Cells', 'Lineage_EVL'
]


adata_1.obs['lineages'] = adata_1.obs[lineage_columns].apply(
    lambda row: ', '.join([col for col in lineage_columns if row[col] == 'TRUE']),
    axis=1
)

# Step 2: Create a new column indicating presence of any lineage


sc.pp.filter_cells(adata_1, min_genes = 200)
sc.pp.filter_cells(adata_1, max_genes=5000)
sc.pp.filter_genes(adata_1, min_cells=3)
sc.pp.normalize_total(adata_1, target_sum=1e4)
sc.pp.log1p(adata_1)
sc.pp.highly_variable_genes(adata_1, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata_1, max_value=10)
sc.tl.pca(adata_1, svd_solver='arpack')
sc.pp.neighbors(adata_1)
sc.tl.umap(adata_1)
sc.tl.louvain(adata_1)


# Clean up 'lineages' column by stripping whitespace
adata_1.obs['lineages'] = adata_1.obs['lineages'].str.strip()

lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)


unwanted_lineages = ['Lineage_Spinal_Cord, Lineage_Diencephalon, Lineage_Optic_Cup, Lineage_Midbrain_Neural_Crest, Lineage_Hindbrain_R3, Lineage_Hindbrain_R4_5_6, Lineage_Telencephalon, Lineage_Epidermis, Lineage_Neural_Plate_Border, Lineage_Placode_Adeno._Lens_Trigeminal, Lineage_Placode_Epibranchial_Otic, Lineage_Placode_Olfactory',
                    ' ', 'Lineage_Tailbud, Lineage_Adaxial_Cells, Lineage_Somites', 'Lineage_Tailbud, Lineage_Adaxial_Cells',
                    'Lineage_Tailbud, Lineage_Adaxial_Cells, Lineage_Somites, Lineage_Hematopoeitic_ICM, Lineage_Hematopoeitic_RBI_Pronephros, Lineage_Endoderm_Pharyngeal, Lineage_Endoderm_Pancreatic_Intestinal, Lineage_Heart_Primordium, Lineage_Cephalic_Mesoderm',
                    'Lineage_Endoderm_Pharyngeal, Lineage_Endoderm_Pancreatic_Intestinal',
                    'Lineage_Spinal_Cord, Lineage_Diencephalon, Lineage_Optic_Cup, Lineage_Midbrain_Neural_Crest, Lineage_Hindbrain_R3, Lineage_Hindbrain_R4_5_6, Lineage_Telencephalon',
                    'Lineage_Hematopoeitic_ICM, Lineage_Hematopoeitic_RBI_Pronephros',
                    'Lineage_Spinal_Cord, Lineage_Diencephalon, Lineage_Optic_Cup, Lineage_Midbrain_Neural_Crest, Lineage_Hindbrain_R3, Lineage_Hindbrain_R4_5_6, Lineage_Telencephalon, Lineage_Epidermis, Lineage_Neural_Plate_Border, Lineage_Placode_Adeno._Lens_Trigeminal, Lineage_Placode_Epibranchial_Otic, Lineage_Placode_Olfactory, Lineage_Tailbud, Lineage_Adaxial_Cells, Lineage_Somites, Lineage_Hematopoeitic_ICM, Lineage_Hematopoeitic_RBI_Pronephros, Lineage_Endoderm_Pharyngeal, Lineage_Endoderm_Pancreatic_Intestinal, Lineage_Heart_Primordium, Lineage_Cephalic_Mesoderm, Lineage_Prechordal_Plate, Lineage_Notochord',
                    'Lineage_Epidermis, Lineage_Neural_Plate_Border',
                    'Lineage_Epidermis, Lineage_Neural_Plate_Border, Lineage_Placode_Adeno._Lens_Trigeminal, Lineage_Placode_Epibranchial_Otic, Lineage_Placode_Olfactory',
                    'Lineage_Placode_Adeno._Lens_Trigeminal',
                    'Lineage_Spinal_Cord, Lineage_Diencephalon',
                    'Lineage_Spinal_Cord, Lineage_Diencephalon, Lineage_Optic_Cup, Lineage_Midbrain_Neural_Crest, Lineage_Hindbrain_R3, Lineage_Hindbrain_R4_5_6, Lineage_Telencephalon, Lineage_Epidermis, Lineage_Neural_Plate_Border, Lineage_Placode_Adeno._Lens_Trigeminal, Lineage_Placode_Epibranchial_Otic, Lineage_Placode_Olfactory, Lineage_Tailbud, Lineage_Adaxial_Cells, Lineage_Somites, Lineage_Hematopoeitic_ICM, Lineage_Hematopoeitic_RBI_Pronephros, Lineage_Endoderm_Pharyngeal, Lineage_Endoderm_Pancreatic_Intestinal, Lineage_Heart_Primordium, Lineage_Cephalic_Mesoderm, Lineage_Prechordal_Plate, Lineage_Notochord, Lineage_Primordial_Germ_Cells',
                    'Lineage_Heart_Primordium, Lineage_Cephalic_Mesoderm',
                    'Lineage_Hematopoeitic_ICM, Lineage_Hematopoeitic_RBI_Pronephros, Lineage_Endoderm_Pharyngeal, Lineage_Endoderm_Pancreatic_Intestinal',
                    'Lineage_Prechordal_Plate, Lineage_Notochord',
                    'Lineage_Spinal_Cord, Lineage_Diencephalon, Lineage_Optic_Cup, Lineage_Midbrain_Neural_Crest']


# Create a dictionary to map unwanted lineages to an empty string
replace_dict = {lineage: '' for lineage in unwanted_lineages}

# Remove unwanted lineages from 'lineages' column
adata_1.obs['lineages'] = adata_1.obs['lineages'].replace(replace_dict, regex=True)


# Remove unwanted lineages from 'lineages' column


# Convert any remaining empty or whitespace-only strings to NaN
adata_1.obs['lineages'].replace('', np.nan, inplace=True)
adata_1.obs['lineages'].replace(' ', np.nan, inplace=True)

# Drop rows with NaN values in 'lineages' column
adata_1 = adata_1[~adata_1.obs['lineages'].isna()]

# Verify the changes
lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)

unwanted_lineages_new = ['Lineage_Tailbud,Lineage_Adaxial_Cells,Lineage_Somites',
                    'Lineage_Tailbud,Lineage_Adaxial_Cells',
                    'Lineage_Tailbud,Lineage_Adaxial_Cells,Lineage_Somites,Lineage_Hematopoeitic_ICM,Lineage_Hematopoeitic_RBI_Pronephros,Lineage_Endoderm_Pharyngeal,Lineage_Endoderm_Pancreatic_Intestinal,Lineage_Heart_Primordium,Lineage_Cephalic_Mesoderm',
                    'Lineage_Endoderm_Pharyngeal,Lineage_Endoderm_Pancreatic_Intestinal',
                    'Lineage_Spinal_Cord,Lineage_Diencephalon,Lineage_Optic_Cup,Lineage_Midbrain_Neural_Crest,Lineage_Hindbrain_R3,Lineage_Hindbrain_R4_5_6,Lineage_Telencephalon',
                    'Lineage_Hematopoeitic_ICM,Lineage_Hematopoeitic_RBI_Pronephros',
                    ',Lineage_Tailbud,Lineage_Adaxial_Cells,Lineage_Somites,Lineage_Hematopoeitic_ICM,Lineage_Hematopoeitic_RBI_Pronephros,Lineage_Endoderm_Pharyngeal,Lineage_Endoderm_Pancreatic_Intestinal,Lineage_Heart_Primordium,Lineage_Cephalic_Mesoderm,Lineage_Prechordal_Plate,Lineage_Notochord',
                    'Lineage_Heart_Primordium,Lineage_Cephalic_Mesoderm',
                    'Lineage_Hematopoeitic_ICM,Lineage_Hematopoeitic_RBI_Pronephros,Lineage_Endoderm_Pharyngeal,Lineage_Endoderm_Pancreatic_Intestinal',
                    'Lineage_Prechordal_Plate,Lineage_Notochord',
                    'Lineage_Spinal_Cord,Lineage_Diencephalon,Lineage_Optic_Cup,Lineage_Midbrain_Neural_Crest']


# Create a dictionary to map unwanted lineages to an empty string
replace_dict = {lineage: '' for lineage in unwanted_lineages_new}

# Remove unwanted lineages from 'lineages' column
adata_1.obs['lineages'] = adata_1.obs['lineages'].replace(replace_dict, regex=True)

# Verify the changes
lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)

# Convert any remaining empty or whitespace-only strings to NaN
adata_1.obs['lineages'].replace('', np.nan, inplace=True)
adata_1.obs['lineages'].replace(' ', np.nan, inplace=True)

# Drop rows with NaN values in 'lineages' column
adata_1 = adata_1[~adata_1.obs['lineages'].isna()]

# Verify the changes
lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)

unwanted_lineages = [',,,', ',,,,,', ',,,,,,Lineage_Primordial_Germ_Cells', ',']

# Create a dictionary to map unwanted lineages to an empty string
replace_dict = {lineage: '' for lineage in unwanted_lineages_new}

# Remove unwanted lineages from 'lineages' column
adata_1.obs['lineages'] = adata_1.obs['lineages'].replace(replace_dict, regex=True)

# Convert any remaining empty or whitespace-only strings to NaN
adata_1.obs['lineages'].replace('', np.nan, inplace=True)
adata_1.obs['lineages'].replace(' ', np.nan, inplace=True)

# Drop rows with NaN values in 'lineages' column
adata_1 = adata_1[~adata_1.obs['lineages'].isna()]


# Verify the changes
lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)


# Remove unwanted lineages from 'lineages' column
for lineage in unwanted_lineages:
    adata_1.obs['lineages'] = adata_1.obs['lineages'].str.replace(lineage, '')


# Remove extra commas and leading/trailing whitespace
adata_1.obs['lineages'] = adata_1.obs['lineages'].str.replace(r',+', ',', regex=True)
adata_1.obs['lineages'] = adata_1.obs['lineages'].str.strip(',')

# Remove empty values
adata_1.obs['lineages'] = adata_1.obs['lineages'].replace('', pd.NA)


# Verify the changes
lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)


# Remove 'Lineage_' prefix from all lineages in the 'lineages' column
adata_1.obs['lineages'] = adata_1.obs['lineages'].str.replace('Lineage_', '')

lineage_counts = adata_1.obs['lineages'].value_counts()
print(lineage_counts)


# Generate the UMAP plot
plt.rcParams.update({'font.size': 18})
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['axes.titlesize'] = 18
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['axes.linewidth'] = 0


sc.pl.umap(adata_1, color=['lineages'], size=20, save="UMAP_lineages_all_zebrafish.pdf")
sc.pl.umap(adata_1, color=['Stage'], size=20, save="UMAP_stages_zebrafish.pdf")
sc.pl.umap(adata_1, color=['Lineage_Cephalic_Mesoderm'], size=20, save="UMAP_Lineage_Cephalic_Mesoderm_zebrafish.pdf")
sc.pl.umap(adata_1, color=['Segment'], size=20, save="UMAP_Segment_zebrafish.pdf")
sc.pl.umap(adata_1, color=['louvain'], size=20, save="UMAP_louvain_zebrafish.pdf")


genes_of_interest = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG",
                     "DKK1", "MESP1", "EOMES", "SALL3", "MIXL1",
                     "STAT3", "FGF8", "FGF10", "FGF4", "ESRRB",
                     "TWIST1", "CDX2", "HAND1", "T"]
# Keep only genes present in the dataset
genes_of_interest_present = [gene for gene in genes_of_interest if gene in adata_1.var.index]

genes_of_interest_present


# Define the custom color map
colors = ["#cccccc", "#ff6666", "#ff0000"]
cmap = mcolors.LinearSegmentedColormap.from_list("", colors)


# Plot the UMAP with the custom color map
sc.pl.umap(adata_1, color=genes_of_interest_present, cmap=cmap, size=5, vmin=0, vmax=10, save="umap_genes_of_interest_zebrafish.pdf")


adata_1_genes_of_interest = adata_1[:, adata_1.var_names.isin(genes_of_interest_present)]

cell_counts = adata_1_genes_of_interest.obs.groupby(['Stage']).size()
print(cell_counts)


# Subset adata to only include cells that express 'Epcam'
adata_1_Epcam = adata_1[adata_1[:, 'EPCAM'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Epcam = adata_1_Epcam.obs.groupby(['Stage']).size()
print(cell_counts_Epcam)


# Subset adata to only include cells that express 'Foxa2'
adata_1_Foxa2 = adata_1[adata_1[:, 'FOXA2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Foxa2 = adata_1_Foxa2.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Gata6'
adata_1_Gata6 = adata_1[adata_1[:, 'GATA6'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Gata6 = adata_1_Gata6.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Wnt11'
adata_1_Wnt11 = adata_1[adata_1[:, 'WNT11'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Wnt11 = adata_1_Wnt11.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Nanog'
adata_1_Nanog = adata_1[adata_1[:, 'NANOG'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Nanog = adata_1_Nanog.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Mixl1'
adata_1_Mixl1 = adata_1[adata_1[:, 'MIXL1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Mixl1 = adata_1_Mixl1.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Stat3'
adata_1_Stat3 = adata_1[adata_1[:, 'STAT3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Stat3 = adata_1_Stat3.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Fgf4'
adata_1_Fgf4 = adata_1[adata_1[:, 'FGF4'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Fgf4 = adata_1_Fgf4.obs.groupby(['Stage']).size()

# Subset adata to only include cells that express 'Esrrb'
adata_1_Esrrb = adata_1[adata_1[:, 'ESRRB'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_Esrrb = adata_1_Esrrb.obs.groupby(['Stage']).size()


# Combine all cell counts into one DataFrame
stages = pd.concat([cell_counts_Epcam, cell_counts_Foxa2, cell_counts_Gata6, cell_counts_Wnt11, cell_counts_Nanog, cell_counts_Mixl1,
                   cell_counts_Stat3, cell_counts_Fgf4, cell_counts_Esrrb], axis=1)
stages.columns = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG",
                     "MIXL1", "STAT3", "FGF4", "ESRRB"]

# Melt the DataFrame to long format for easier plotting
stages_melted = stages.reset_index().melt(id_vars='Stage', var_name='gene', value_name='cell_count')


# Define the gene groups
genes_of_interest_present = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG", "MIXL1", "STAT3", "FGF4", "ESRRB"]

# Define stage order manually
stage_order = ['03.3-HIGH', '03.8-OBLONG', '04.3-DOME', '04.8-30%', '05.3-50%', '06.0-SHIELD', '07.0-60%', '08.0-75%', '09.0-90%', '10.0-BUD', '11.0-3-Somite', '12.0-6-Somite']


# Convert 'stage' to an ordered categorical type
stages_melted['Stage'] = pd.Categorical(stages_melted['Stage'], categories=stage_order, ordered=True)

# Convert 'stage' to numerical values for smoothing
stages_melted['Stage_num'] = stages_melted['Stage'].str.extract(r'(\d+.\d+)').astype(float)

# Sort data for rolling mean
stages_melted.sort_values(by=['gene', 'Stage_num'], inplace=True)

# Compute rolling mean with a window of 3
stages_melted['cell_count_smooth'] = stages_melted.groupby('gene')['cell_count'].transform(lambda x: x.rolling(3, 1).mean())

# Create a single subplot
fig, ax = plt.subplots(figsize=(8, 6))

# Filter the data for the genes of interest
stages_melted_interest = stages_melted[stages_melted['gene'].isin(genes_of_interest_present)]

# Create the line plot for the genes of interest with smoothed cell counts
sns.lineplot(x='Stage', y='cell_count_smooth', hue='gene', data=stages_melted_interest, ax=ax)

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
plt.savefig('Zebrafish_Cell_count_plot_Stages.pdf', format='pdf')

# Show the plot
plt.show()


## Cell counts lineages
cell_counts_lineages = adata_1_genes_of_interest.obs.groupby(['lineages']).size()
print(cell_counts_lineages)


# Subset adata to only include cells that express 'Epcam'
adata_1_lineage_Epcam = adata_1[adata_1[:, 'EPCAM'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Epcam = adata_1_lineage_Epcam.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Foxa2'
adata_1_lineage_Foxa2 = adata_1[adata_1[:, 'FOXA2'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Foxa2 = adata_1_lineage_Foxa2.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Gata6'
adata_1_lineage_Gata6 = adata_1[adata_1[:, 'GATA6'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Gata6 = adata_1_lineage_Gata6.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Wnt11'
adata_1_lineage_Wnt11 = adata_1[adata_1[:, 'WNT11'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Wnt11 = adata_1_lineage_Wnt11.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Nanog'
adata_1_lineage_Nanog = adata_1[adata_1[:, 'NANOG'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Nanog = adata_1_lineage_Nanog.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Mixl1'
adata_1_lineage_Mixl1 = adata_1[adata_1[:, 'MIXL1'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Mixl1 = adata_1_lineage_Mixl1.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Stat3'
adata_1_lineage_Stat3 = adata_1[adata_1[:, 'STAT3'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Stat3 = adata_1_lineage_Stat3.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Fgf4'
adata_1_lineage_Fgf4 = adata_1[adata_1[:, 'FGF4'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Fgf4 = adata_1_lineage_Fgf4.obs.groupby(['lineages']).size()

# Subset adata to only include cells that express 'Esrrb'
adata_1_lineage_Esrrb = adata_1[adata_1[:, 'ESRRB'].X > 0]
# Count the number of cells in each stage for each cluster
cell_counts_lineage_Esrrb = adata_1_lineage_Esrrb.obs.groupby(['lineages']).size()


# Combine all cell counts into one DataFrame
celltypes = pd.concat([cell_counts_lineage_Epcam, cell_counts_lineage_Foxa2, cell_counts_lineage_Gata6,
                       cell_counts_lineage_Wnt11, cell_counts_lineage_Nanog,
                       cell_counts_lineage_Mixl1, cell_counts_lineage_Stat3,
                       cell_counts_lineage_Fgf4, cell_counts_lineage_Esrrb], axis=1)
celltypes.columns = ["EPCAM", "FOXA2", "GATA6", "WNT11", "NANOG", "MIXL1",
                     "STAT3", "FGF4", "ESRRB"]

# Melt the DataFrame to long format for easier plotting
celltypes_melted = celltypes.reset_index().melt(id_vars='lineages', var_name='gene', value_name='cell_count')


# Define size factor for your dots. You might have to adjust this value to get your desired output.
size_factor = 4.5

# Create a scatter plot
plt.figure(figsize=(30, 30))
scatter_plot = sns.scatterplot(x='gene', y='lineages', size='cell_count',
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
plt.savefig('Zebrafish_dotplot_lineages_cell_counts.pdf', format='pdf', bbox_inches='tight')

# Show the plot
plt.show()


Epcam_Foxa2 = ["EPCAM", "FOXA2"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2 = [adata_1[:, gene].X > 0 for gene in Epcam_Foxa2]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2 = np.all(masks_Epcam_Foxa2, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Foxa2'] = coexpressed_mask_Epcam_Foxa2.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Foxa2', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Foxa2.pdf")


Epcam_Gata6 = ["EPCAM", "GATA6"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Gata6 = [adata_1[:, gene].X > 0 for gene in Epcam_Gata6]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Gata6 = np.all(masks_Epcam_Gata6, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Gata6'] = coexpressed_mask_Epcam_Gata6.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Gata6', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Gata6.pdf")


Epcam_Wnt11 = ["EPCAM", "WNT11"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Wnt11 = [adata_1[:, gene].X > 0 for gene in Epcam_Wnt11]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Wnt11 = np.all(masks_Epcam_Wnt11, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Wnt11'] = coexpressed_mask_Epcam_Wnt11.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Wnt11', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Wnt11.pdf")


Foxa2_Gata6 = ["FOXA2", "GATA6"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Foxa2_Gata6 = [adata_1[:, gene].X > 0 for gene in Foxa2_Gata6]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Foxa2_Gata6 = np.all(masks_Foxa2_Gata6, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Foxa2_Gata6'] = coexpressed_mask_Foxa2_Gata6.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Foxa2_Gata6', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Foxa2_Gata6.pdf")


Epcam_Foxa2_Gata6 = ["EPCAM", "FOXA2", "GATA6"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2_Gata6 = [adata_1[:, gene].X > 0 for gene in Epcam_Foxa2_Gata6]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2_Gata6 = np.all(masks_Epcam_Foxa2_Gata6, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Foxa2_Gata6'] = coexpressed_mask_Epcam_Foxa2_Gata6.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Foxa2_Gata6', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Foxa2_Gata6.pdf")


Epcam_Foxa2_Gata6_Wnt11 = ["EPCAM", "FOXA2", "GATA6", "WNT11"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Foxa2_Gata6_Wnt11 = [adata_1[:, gene].X > 0 for gene in Epcam_Foxa2_Gata6_Wnt11]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Foxa2_Gata6_Wnt11 = np.all(masks_Epcam_Foxa2_Gata6_Wnt11, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Foxa2_Gata6_Wnt11'] = coexpressed_mask_Epcam_Foxa2_Gata6_Wnt11.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Foxa2_Gata6_Wnt11', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Foxa2_Gata6_Wnt11.pdf")


Epcam_Stat3 = ["EPCAM", "STAT3"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Stat3 = [adata_1[:, gene].X > 0 for gene in Epcam_Stat3]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Stat3 = np.all(masks_Epcam_Stat3, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Stat3'] = coexpressed_mask_Epcam_Stat3.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Stat3', cmap=cmap, size=5,
           save="_Zebrafish_coexpressed_Epcam_Stat3.pdf")


Epcam_Gata6_Wnt11 = ["EPCAM", "GATA6", "WNT11"]

# Create a boolean mask for each gene indicating whether it is expressed in each cell
masks_Epcam_Gata6_Wnt11 = [adata_1[:, gene].X > 0 for gene in Epcam_Gata6_Wnt11]

# Combine the masks using logical AND to get a mask indicating whether all genes are expressed
coexpressed_mask_Epcam_Gata6_Wnt11 = np.all(masks_Epcam_Gata6_Wnt11, axis=0)

# Convert the boolean mask to integer
adata_1.obs['coexpressed_Epcam_Gata6_Wnt11'] = coexpressed_mask_Epcam_Gata6_Wnt11.astype(int)

# Plot the UMAP highlighting cells where all signature genes are co-expressed
sc.pl.umap(adata_1, color='coexpressed_Epcam_Gata6_Wnt11', cmap=cmap, size=5, save="_Zebrafish_coexpressed_Epcam_Gata6_Wnt11.pdf")


# Subset adata to only include cells that coexpress "Epcam" and "Foxa2"
adata_celltype_Epcam_Foxa2 = adata_1[np.logical_and(adata_1[:, "EPCAM"].X > 0,
                                                           adata_1[:, "FOXA2"].X > 0)]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2 = adata_celltype_Epcam_Foxa2.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Gata6"
adata_celltype_Epcam_Gata6 = adata_1[np.logical_and(adata_1[:, "EPCAM"].X > 0, adata_1[:, "GATA6"].X > 0)]
# Count the number of cells in each stage
cell_counts_Epcam_Gata6 = adata_celltype_Epcam_Gata6.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Epcam" and "Wnt11"
adata_celltype_Epcam_Wnt11 = adata_1[np.logical_and(adata_1[:, "EPCAM"].X > 0, adata_1[:, "WNT11"].X > 0)]
# Count the number of cells in each stage
cell_counts_Epcam_Wnt11 = adata_celltype_Epcam_Wnt11.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Foxa2" and "Gata6"
adata_celltype_Foxa2_Gata6 = adata_1[np.logical_and(adata_1[:, "FOXA2"].X > 0, adata_1[:, "GATA6"].X > 0)]
# Count the number of cells in each stage
cell_counts_Foxa2_Gata6 = adata_celltype_Foxa2_Gata6.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Foxa2", and "Gata6"
adata_celltype_Epcam_Foxa2_Gata6 = adata_1[np.logical_and.reduce((adata_1[:, "EPCAM"].X > 0,
                                                                  adata_1[:, "FOXA2"].X > 0,
                                                                  adata_1[:, "GATA6"].X > 0))]
# Count the number of cells in each stage
cell_counts_Epcam_Foxa2_Gata6 = adata_celltype_Epcam_Foxa2_Gata6.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Foxa2", "Gata6", and "Wnt11"
adata_celltype_Epcam_Foxa2_Gata6_Wnt11 = adata_1[
    np.logical_and.reduce(
        (
            adata_1[:, "EPCAM"].X > 0,
            adata_1[:, "FOXA2"].X > 0,
           adata_1[:, "GATA6"].X > 0,
            adata_1[:, "WNT11"].X > 0,
        )
    )
]

# Count the number of cells in each stage
cell_counts_Epcam_Foxa2_Gata6_Wnt11 = adata_celltype_Epcam_Foxa2_Gata6_Wnt11.obs.groupby(["Stage"]).size()


# Subset adata to only include cells that coexpress "Epcam" and "Stat3"
adata_celltype_Epcam_Stat3 = adata_1[np.logical_and.reduce((adata_1[:, "EPCAM"].X > 0, adata_1[:, "STAT3"].X > 0))]

# Count the number of cells in each stage
cell_counts_Epcam_Stat3 = adata_celltype_Epcam_Stat3.obs.groupby(['Stage']).size()


# Subset adata to only include cells that coexpress "Epcam", "Gata6", and "Wnt11"
adata_celltype_Epcam_Gata6_Wnt11 = adata_1[np.logical_and.reduce((adata_1[:, "EPCAM"].X > 0,
                                                                  adata_1[:, "GATA6"].X > 0,
                                                                  adata_1[:, "WNT11"].X > 0))]

# Count the number of cells in each stage
cell_counts_Epcam_Gata6_Wnt11 = adata_celltype_Epcam_Gata6_Wnt11.obs.groupby(['Stage']).size()


## Coexpression data merging
merged_data = pd.concat([cell_counts_Epcam_Foxa2,
                         cell_counts_Epcam_Gata6, cell_counts_Epcam_Wnt11, cell_counts_Foxa2_Gata6,
                         cell_counts_Epcam_Foxa2_Gata6, cell_counts_Epcam_Foxa2_Gata6_Wnt11,
                        cell_counts_Epcam_Stat3, cell_counts_Epcam_Gata6_Wnt11], axis=1)


merged_data.columns = ['Epcam_Foxa2', 'Epcam_Gata6', 'Epcam_Wnt11', 'Foxa2_Gata6', 'Epcam_Foxa2_Gata6', 'Epcam_Foxa2_Gata6_Wnt11',
                       'Epcam_Stat3', 'Epcam_Gata6_Wnt11']


merged_data.fillna(0, inplace=True)

merged_data.to_csv('Zebrafish_coexpression_merged_data_stages.csv', index=False)


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
plt.savefig('Zebrafish_Lineplot_Coexpression_singlecells_count.pdf', format='pdf')

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
plt.savefig('Zebrafish_coexpression_new_heatmap.pdf', format='pdf')

# Show the plot
plt.show()
