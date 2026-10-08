# Data

No data are stored in this repository. Download the datasets below and place the files in the
sub-folders indicated (folder names are those expected by the scripts). All paths are relative
to the repository root.

| Folder | Dataset | Source |
|---|---|---|
| `data/bulk_RNAseq/` | 3′-RNA-seq, WT and *Epcam*<sup>−/−</sup> mESC (#56, #114), D0–D10 | GEO [GSE293121](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE293121) (this study) |
| `data/Mohammed2017_GSE100597/` | Mouse embryos E3.5–E6.75 (Smart-seq2) | GEO [GSE100597](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100597) |
| `data/PijuanSala2019_atlas/` | Mouse gastrulation atlas E6.5–E8.5 (10x) | ArrayExpress [E-MTAB-6967](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6967); processed atlas files from the [MarioniLab/EmbryoTimecourse2018](https://github.com/MarioniLab/EmbryoTimecourse2018) repository |
| `data/Dong2018_GSE87038/` | Mouse organogenesis E9.5–E11.5 | GEO [GSE87038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87038) |
| `data/Xu2023_GSE157329/` | Human embryos CS12–CS16 (10x) | GEO [GSE157329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157329) |
| `data/Farrell2018_zebrafish_URD/` | Zebrafish embryogenesis, Drop-seq URD atlas | GEO [GSE106587](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106587) |

The single-cell RNA-seq data of WT mESC embryoid bodies generated in this study are available
under GEO [GSE318465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318465).

---

## Expected files

### `data/bulk_RNAseq/`

| File | Description |
|---|---|
| `featurecounts_mESC_EpCAM_bulkRNAseq.Rdata` | `.Rdata` file containing a single object: the output list of `Rsubread::featureCounts()` (element `$counts` = Ensembl gene IDs × samples). |

The scripts derive the sample groups from the column names, which follow the pattern
`Sample.<GROUP>.<DAY>.<REPLICATE>...Aligned.sortedByCoord.out.bam`, e.g.
`Sample.WT.D3.2.1Aligned.sortedByCoord.out.bam`, with `<GROUP>` ∈ {`WT`, `N56`, `N114`},
`<DAY>` ∈ {`D0`, `D3`, `D7`, `D10`} and `<REPLICATE>` ∈ {1–4}.
To start from a count matrix (genes × samples) instead, wrap it in the same structure:

```r
fc <- list(counts = as.matrix(count_matrix))   # colnames as above
save(fc, file = "data/bulk_RNAseq/featurecounts_mESC_EpCAM_bulkRNAseq.Rdata")
```

### `data/Mohammed2017_GSE100597/`

| File | Description |
|---|---|
| `GSE100597_count_table_QC_filtered.txt` | Supplementary file of GSE100597 (gunzip). Genes × cells; gene IDs `Symbol-chr-start`; cell names start with the stage (`E3.5_…`, `E4.5_…`, `E5.5_…`, `E6.5_…`, `E6.75_…`). |

### `data/PijuanSala2019_atlas/`

| File | Description |
|---|---|
| `raw_counts.mtx`, `genes.tsv`, `barcodes.tsv`, `meta.tab` | Processed atlas files (`atlas_data.tar.gz`, see the MarioniLab/EmbryoTimecourse2018 repository; also available through the Bioconductor package `MouseGastrulationData`). |
| `adata_2.txt` | Tab-separated expression table (cells/genes) exported from the atlas object; input of `03_atlas_subset_Seurat_marker_genes.R`. |

### `data/Dong2018_GSE87038/`

| File | Description |
|---|---|
| `GSE87038_Mouse_Organogenesis_UMI_counts_matrix.txt` | Genes × cells UMI count matrix (supplementary file of GSE87038). |
| `heart.txt` | Subset of the matrix above containing the heart-derived cells (first column = gene). |

### `data/Xu2023_GSE157329/`

| File | Description |
|---|---|
| `GSE157329_raw_counts.mtx` | Sparse count matrix (genes × cells). |
| `GSE157329_gene_annotate.txt` | Gene annotation (`gene_id`, `gene_short_name`). |
| `GSE157329_cell_annotate_1.txt` | `GSE157329_cell_annotate.txt` with the `barcode` column moved to the first position. |
| `barcode_new_1.txt` | One-column file (header `barcode`) with the cell barcodes in the column order of the count matrix. |

### `data/Farrell2018_zebrafish_URD/`

| File | Description |
|---|---|
| `URD_Dropseq_Expression_Log2TPM.txt` | Genes × cells log2(TPM) expression (first column `GENE`). |
| `meta.txt` | Cell metadata (`NAME`, `Stage`, `Segment`, `Lineage_*` flags). |

---

## References for the public datasets

- Mohammed H, *et al.* Single-cell landscape of transcriptional heterogeneity and cell fate decisions during mouse early gastrulation. *Cell Rep* 2017;20:1215–1228.
- Pijuan-Sala B, *et al.* A single-cell molecular map of mouse gastrulation and early organogenesis. *Nature* 2019;566:490–495.
- Dong J, *et al.* Single-cell RNA-seq analysis unveils a prevalent epithelial/mesenchymal hybrid state during mouse organogenesis. *Genome Biol* 2018;19:31.
- Xu Y, *et al.* A single-cell transcriptome atlas profiles early organogenesis in human embryos. *Nat Cell Biol* 2023;25:604–615.
- Farrell JA, *et al.* Single-cell reconstruction of developmental trajectories during zebrafish embryogenesis. *Science* 2018;360:eaar3131.
