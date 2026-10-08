# =============================================================================
# Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
#           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
# Module  : 03 | Public scRNA-seq: mouse gastrulation atlas E6.5-E8.5
#           (Pijuan-Sala et al., Nature 2019; ArrayExpress E-MTAB-6967) and
#           mouse organogenesis E9.5-E11.5 (Dong et al., Genome Biol 2018; GEO: GSE87038)
# Script  : 01_atlas_preprocessing_stage_pseudobulk.py
# Purpose : Scanpy preprocessing of the Pijuan-Sala et al. (2019) gastrulation atlas: QC,
#           normalisation, HVG, PCA, Louvain, UMAP; annotation with sample/stage metadata and
#           stage-level pseudobulk matrix (summed counts per stage) with UMAP/Louvain.
# Input   : data/PijuanSala2019_atlas/{raw_counts.mtx, genes.tsv, barcodes.tsv, meta.tab}
# Output  : results/03_PijuanSala2019/{gene_matrix.csv, stage_matrix.csv}
# Notes   : Converted from the original Jupyter notebook (jupytext 'percent' format: one
#           # %% per cell). gene_matrix.csv densifies the full atlas and requires a lot of RAM.
# =============================================================================

# %%
import os
from pathlib import Path

# ---- Paths --------------------------------------------------------------------
# Run from the repository root. Inputs in data/PijuanSala2019_atlas/ (see
# data/README.md); outputs are written to results/03_PijuanSala2019/.
DATA_DIR = Path("data/PijuanSala2019_atlas").resolve()
OUT_DIR = Path("results/03_PijuanSala2019").resolve()
OUT_DIR.mkdir(parents=True, exist_ok=True)
os.chdir(OUT_DIR)

# %%
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.io import mmread
import anndata
import anndata as ad
import scipy.io

# %%
from sklearn.preprocessing import StandardScaler
import umap

# %%
adata = sc.read_mtx(str(DATA_DIR / "raw_counts.mtx"))

# %%
adata = adata.T

# %%
# Read the gene information file
genes_df = pd.read_csv(str(DATA_DIR / "genes.tsv"), header=None, sep='\t', index_col=0, names=['gene_ids', 'gene_symbols'])

# %%
genes_df

# %%
# Set the index for adata.var and map gene_symbols using gene_ids
adata.var.index = genes_df.index
adata.var['gene_symbols'] = adata.var.index.map(genes_df['gene_symbols'])

# %%
barcodes = pd.read_csv(str(DATA_DIR / "barcodes.tsv"), header=None, sep='\t', names=['cell_barcode'])

# %%
adata.obs.index = barcodes['cell_barcode']

# %%
dense_matrix = adata.X.toarray()

# %%
dense_matrix

# %%
gene_matrix = pd.DataFrame(dense_matrix, index=adata.obs.index, columns=adata.var['gene_symbols'])

# %%
gene_matrix

# %%
gene_matrix.to_csv("gene_matrix.csv")

# %%
print("adata shape:", adata.shape)
print("adata.var shape:", adata.var.shape)

# %%
print(adata.var.head())

# %%
num_genes_per_cell = (gene_matrix > 0).sum(axis=1)

# %%
UMIs_per_cell = gene_matrix.sum(axis=1)

# %%
adata.obs['n_genes'] = num_genes_per_cell
adata.obs['n_counts'] = UMIs_per_cell

# %%
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sns.histplot(adata.obs['n_genes'], bins=50, ax=axes[0], kde=False)
axes[0].set_xlabel('Number of genes')
axes[0].set_ylabel('Number of cells')
axes[0].set_title('Number of genes per cell')

sns.histplot(adata.obs['n_counts'], bins=50, ax=axes[1], kde=False)
axes[1].set_xlabel('Number of counts')
axes[1].set_ylabel('Number of cells')
axes[1].set_title('Number of counts per cell')

plt.tight_layout()
plt.show()

# %%
adata

# %%
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# %%
sc.pp.normalize_total(adata, target_sum=1e4)

# %%
sc.pp.log1p(adata)

# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# %%
sc.tl.pca(adata, svd_solver='arpack')

# %%
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# %%
sc.tl.louvain(adata)

# %%
sc.tl.umap(adata)
sc.pl.umap(adata, color=['louvain'], legend_loc='on data', title='', frameon=False, palette='tab20')

# %%
metadata = pd.read_csv(str(DATA_DIR / "meta.tab"), sep='\t')

# %%
print(metadata.head())

# %%
print(metadata.columns)

# %%
adata.obs = adata.obs.merge(metadata[['cell', 'sample']], left_index=True, right_on='cell', how='left')

# %%
sample_matrix = gene_matrix.groupby(adata.obs['sample']).sum()

# %%
sample_matrix

# %%
print(adata.obs.head())

# %%
print("adata.obs index (cell barcodes):")
print(adata.obs.index[:5])

print("\nmetadata cell barcodes:")
print(metadata['cell'].head())

# %%
print("Unique sample identifiers:")
print(adata.obs['sample'].unique())

# %%
metadata['cell'] = metadata['cell'].str.replace('cell_', '').astype(adata.obs.index.dtype)

# %%
adata.obs = adata.obs.merge(metadata[['cell', 'stage']], left_index=True, right_on='cell', how='left')

# %%
stage_matrix = gene_matrix.groupby(adata.obs['stage']).sum()

# %%
stage_matrix

# %%
print("Unique stages:")
print(adata.obs['stage'].unique())

# %%
gene_matrix['stage'] = adata.obs['stage'].values

# %%
stage_matrix = gene_matrix.groupby('stage').sum()

# %%
print(stage_matrix.shape)
print(stage_matrix.head())

# %%
scaler = StandardScaler()
stage_matrix_scaled = scaler.fit_transform(stage_matrix)

# %%
reducer = umap.UMAP(random_state=42)
stage_embedding = reducer.fit_transform(stage_matrix_scaled)

# %%
stage_to_num = {stage: i for i, stage in enumerate(adata.obs['stage'].unique())}
adata.obs['stage_num'] = adata.obs['stage'].map(stage_to_num)

# %%
adata.obs.reset_index(drop=True, inplace=True)

# %%
stages = ['E6.5', 'E6.75', 'E7.0', 'E7.25', 'E7.5', 'E7.75', 'E8.0', 'E8.25', 'E8.5', 'mixed_gastrulation']

# %%
stage_matrix

# %%
stage_matrix

# %%
adata.obs['stage'] = pd.Categorical(adata.obs['stage'])

# %%
stage_matrix.to_csv("stage_matrix.csv")

# %%
print(stage_matrix.head())

# %%
# Create an AnnData object
adata = sc.AnnData(stage_matrix)

# %%
# Remove NaN values in observation metadata
adata.obs = adata.obs.dropna()

# Remove NaN values in variable metadata
adata.var = adata.var.dropna()

# %%
# Preprocess the data
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata)

# %%
# Compute UMAP
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# %%
# Compute Louvain clusters
sc.tl.louvain(adata)

# %%
# Plot the UMAP with Louvain clusters
sc.pl.umap(adata, color=['louvain'])
