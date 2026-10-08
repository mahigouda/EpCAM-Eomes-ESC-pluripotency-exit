## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 01 | Bulk 3'-RNA-seq of WT and Epcam-/- mESC (clones #56, #114), D0-D10 (GEO: GSE293121)
## Script  : 02_DESeq2_DEGs_lfc1_target_genes.R
## Purpose : Primary differential expression analysis (DESeq2 + IHW, lfcThreshold = 1).
##           Pairwise KO-vs-WT contrasts per day, time-course designs, pooled KO (#56/#114)
##           vs WT per day, volcano plots, top-DEG heatmaps, GO over-representation and
##           expression of selected target genes (e.g. Eomes, Foxa2, Gata6, Hand1, Cd86).
## Input   : data/bulk_RNAseq/featurecounts_mESC_EpCAM_bulkRNAseq.Rdata
## Output  : results/01_bulk_RNAseq/  (Threshold_1_*.xlsx DE tables, target-gene tables)
## Notes   : Sample groups: WT, N56 (= KO clone #56), N114 (= KO clone #114) x D0/D3/D7/D10,
##           n = 4 per group. Run interactively in RStudio for plots.
## =============================================================================

## ---- Packages -----------------------------------------------------------------

library(Rsubread)
library(DESeq2)
library(EnsDb.Mmusculus.v79)
library(EnsDb.Hsapiens.v86)
library(pcaExplorer)
library(IHW)
library(pheatmap)
library(gplots)
library(EnhancedVolcano)
library(WriteXLS)
library(DT)
library(readxl)
library(tidyverse)
library(gdata)
library(knitr)
library(gtools)
library(openxlsx)
library(BiocParallel)
library(purrr)
library(grid)
library(MatrixGenerics)
library(matrixStats)
library(ShortRead)
library(GenomicAlignments)
library(Biobase)
library(GenomicRanges)
library(Biostrings)
library(GenomeInfoDb)
library(XVector)
library(IRanges)
library(S4Vectors)
library(BiocGenerics)
library(dplyr)
library(ggplot2)
library(splines)
library(VGAM)
library(scater)
library(Rsamtools)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggrepel)
library(MAGeCKFlute)
library(pathview)
library(tidyr)
library(reshape2)
library(depmap)
library(ExperimentHub)

## ---- Paths ------------------------------------------------------------------
## Run this script from the repository root. Input files are expected in
## data/bulk_RNAseq/ (see data/README.md); all outputs are written to
## results/01_bulk_RNAseq/.
data_dir <- file.path(getwd(), "data", "bulk_RNAseq")
out_dir  <- file.path(getwd(), "results", "01_bulk_RNAseq")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

## ---- Import featureCounts, annotate genes, build DESeq2 object ------------

## featureCounts output object (list with element $counts: Ensembl gene IDs x samples)
fc.mouse <- mget(load(file.path(data_dir, "featurecounts_mESC_EpCAM_bulkRNAseq.Rdata")))[[1]]
fc.mm <- fc.mouse$counts

smps <- colnames(fc.mm)
smps <- gsub("Sample.","",smps)
smps <- gsub("Aligned.sortedByCoord.out.bam","",smps)
smps <- gsub(".1$","",smps)
cs.mm <- colSums(fc.mm)
names(cs.mm) <- smps

smps.r <- gsub("\\.1","",smps)
smps.r <- gsub("\\.2","",smps.r)
smps.r <- gsub("\\.3","",smps.r)
smps.r <- gsub("\\.4","",smps.r)

table(smps.r)
smps.r

cs.df <- data.frame(names=c(smps),mreads = c( cs.mm), samplenames = smps.r)

ggplot(cs.df, aes(x=names, y=mreads,colour = samplenames))+ geom_point(size = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylim(c(0,10^7))

keep.mm <- rowSums(fc.mm) > 100

table(keep.mm)
keep.mm

fc.mm.filt <- fc.mm[keep.mm,]
sms <- colnames(fc.mm.filt)
sms <- gsub("Sample.","",sms)
sms <- gsub(".1Aligned.sortedByCoord.out.bam","",sms)
colnames(fc.mm.filt) <- sms
ens.cts.mm <- row.names(fc.mm.filt)

anno.result.mm <- AnnotationDbi::select(EnsDb.Mmusculus.v79, keys=ens.cts.mm, columns=c("GENEID","SYMBOL","GENENAME","ENTREZID"),keytype="GENEID")
fc.mm.filt <- fc.mm.filt[anno.result.mm$GENEID,]

row.names(fc.mm.filt) <- anno.result.mm$SYMBOL

hd.pc.pca.mm.dat <- data.frame(gene=row.names(fc.mm.filt),fc.mm.filt)
hd.pc.pca.mm.coll <- aggregate(. ~ gene, data = hd.pc.pca.mm.dat, sum)
row.names(hd.pc.pca.mm.coll) <- hd.pc.pca.mm.coll$gene
hd.pc.pca.mm.coll <- hd.pc.pca.mm.coll[-1,-1]
em.mm <- hd.pc.pca.mm.coll
treatment <- smps.r
treatment <- factor(treatment)

fakeSampleTable.mm <- data.frame(sampleNames = sms, sampleFiles = sms, treatment = treatment)
ddsHTSeq.mm <- DESeqDataSetFromMatrix(countData = em.mm, colData  = fakeSampleTable.mm,  design= ~treatment)
vst.mm <- vst(ddsHTSeq.mm)
pcaplot(vst.mm,intgroup = c("treatment"), ntop = 10000,
        pcX = 1, pcY = 2, title = "PCA RNAseq mouse", text_labels = T, ellipse.prob = .5)

pcaplot(vst.mm,intgroup = c("treatment"), ntop = 10000,)

pcaobj.mm <- prcomp(t(assay(vst.mm)))

hi_loadings(pcaobj.mm,topN = 40, title = "Top/bottom loadings mouse ")
hi_loadings(pcaobj.mm,topN = 40, title = "Top/bottom loadings mouse topN 40 ")
hi_loadings(pcaobj.mm,topN = 40, title = "Top/bottom loadings mouse topN 100 ")
hi_loadings(pcaobj.mm,topN = 100, title = "Top/bottom loadings mouse topN 100 ")
a.vst.mm <- assay(vst.mm)
var.mm <- apply(a.vst.mm, MARGIN = 1, FUN = var)
o.var.mm <- order(var.mm, decreasing = T)
pheatmap(a.vst.mm[o.var.mm,][1:50,], main = "mouse")

## ---- Pairwise DE: each Epcam-/- clone vs WT at D0, D3, D7, D10 (c1-c8) -----

ddsHTSeq.mm.c1 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D0","WT.D0")]

ddsHTSeq.mm.c1$treatment <- factor(as.character(ddsHTSeq.mm.c1$treatment))
ddsHTSeq.mm.c1$treatment <- relevel(ddsHTSeq.mm.c1$treatment, ref = "WT.D0")
dds.c1 <- DESeq(ddsHTSeq.mm.c1)

res.c1 <- results(dds.c1, filterFun=ihw, lfcThreshold = 1)

summary(res.c1)

res.c1.plot <-  res.c1[!is.na(res.c1$padj),]
v.c1 <- data.frame(res.c1.plot[,c(2,6)])
colnames(v.c1)[2] <- "FDR"
ev.c1 <- EnhancedVolcano(v.c1,
                         lab = rownames(v.c1),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N114 D0 Vs WT D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))

ev.c1
res.sig.c1 <- res.c1[!is.na(res.c1$padj),]
res.sig.c1 <- res.sig.c1[abs(res.sig.c1$log2FoldChange) > 1,]
res.sig.c1 <- res.sig.c1[res.sig.c1$padj < 0.1,]
res.sig.c1 <- res.sig.c1[order(res.sig.c1$pvalue),]
res.c1.df <- data.frame(res.c1)
res.sig.c1.df <- data.frame(res.sig.c1)

WriteXLS(c("res.c1.df","res.sig.c1.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N114.D0_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)
ddsHTSeq.mm.c1 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","WT.D3")]

ddsHTSeq.mm.c1$treatment <- factor(as.character(ddsHTSeq.mm.c1$treatment))
ddsHTSeq.mm.c1$treatment <- relevel(ddsHTSeq.mm.c1$treatment, ref = "WT.D3")
ddsHTSeq.mm.c2 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","WT.D3")]

ddsHTSeq.mm.c2$treatment <- factor(as.character(ddsHTSeq.mm.c2$treatment))
ddsHTSeq.mm.c2$treatment <- relevel(ddsHTSeq.mm.c1$treatment, ref = "WT.D3")
dds.c2 <- DESeq(ddsHTSeq.mm.c2)

res.c2 <- results(dds.c2, filterFun=ihw, lfcThreshold = 1)

summary(res.c2)

res.c2.plot <-  res.c2[!is.na(res.c2$padj),]
v.c2 <- data.frame(res.c2.plot[,c(2,6)])
colnames(v.c2)[2] <- "FDR"
ev.c2 <- EnhancedVolcano(v.c2,
                         lab = rownames(v.c2),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N114 D3 Vs WT D3',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c2
res.sig.c2 <- res.c2[!is.na(res.c2$padj),]
res.sig.c2 <- res.sig.c2[abs(res.sig.c2$log2FoldChange) > 1,]
res.sig.c2 <- res.sig.c2[res.sig.c2$padj < 0.1,]
res.sig.c2 <- res.sig.c2[order(res.sig.c2$pvalue),]
res.c2.df <- data.frame(res.c2)
res.sig.c2.df <- data.frame(res.sig.c2)

WriteXLS(c("res.c2.df","res.sig.c2.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N114.D3_vs_WT.D3.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c3 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D7","WT.D7")]
ddsHTSeq.mm.c3$treatment <- factor(as.character(ddsHTSeq.mm.c3$treatment))
ddsHTSeq.mm.c3$treatment <- relevel(ddsHTSeq.mm.c3$treatment, ref = "WT.D7")
dds.c3 <- DESeq(ddsHTSeq.mm.c3)

res.c3 <- results(dds.c3, filterFun=ihw, lfcThreshold = 1)

summary(res.c3)

res.c3.plot <-  res.c3[!is.na(res.c3$padj),]

v.c3 <- data.frame(res.c3.plot[,c(2,6)])
colnames(v.c3)[2] <- "FDR"
ev.c3 <- EnhancedVolcano(v.c3,
                         lab = rownames(v.c3),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N114 D7 Vs WT D7',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c3

res.sig.c3 <- res.c3[!is.na(res.c3$padj),]
res.sig.c3 <- res.sig.c3[abs(res.sig.c3$log2FoldChange) > 1,]
res.sig.c3 <- res.sig.c3[res.sig.c3$padj < 0.1,]
res.sig.c3 <- res.sig.c3[order(res.sig.c3$pvalue),]
res.c3.df <- data.frame(res.c3)
res.sig.c3.df <- data.frame(res.sig.c3)
WriteXLS(c("res.c3.df","res.sig.c3.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N114.D7_vs_WT.D7.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c4 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D10","WT.D10")]
ddsHTSeq.mm.c4$treatment <- factor(as.character(ddsHTSeq.mm.c4$treatment))
ddsHTSeq.mm.c4$treatment <- relevel(ddsHTSeq.mm.c4$treatment, ref = "WT.D10")
dds.c4 <- DESeq(ddsHTSeq.mm.c4)

res.c4 <- results(dds.c4, filterFun=ihw, lfcThreshold = 1)

summary(res.c4)

res.c4.plot <-  res.c4[!is.na(res.c4$padj),]

v.c4 <- data.frame(res.c4.plot[,c(2,6)])
colnames(v.c4)[2] <- "FDR"
ev.c4 <- EnhancedVolcano(v.c4,
                         lab = rownames(v.c4),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N114 D10 Vs WT D10 versus D10',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c4

res.sig.c4 <- res.c4[!is.na(res.c4$padj),]
res.sig.c4 <- res.sig.c4[abs(res.sig.c4$log2FoldChange) > 1,]
res.sig.c4 <- res.sig.c4[res.sig.c4$padj < 0.1,]
res.sig.c4 <- res.sig.c4[order(res.sig.c4$pvalue),]
res.c4.df <- data.frame(res.c4)
res.sig.c4.df <- data.frame(res.sig.c4)

WriteXLS(c("res.c4.df","res.sig.c4.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N114.D10_vs_WT.D10.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c5 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D0","WT.D0")]
ddsHTSeq.mm.c5$treatment <- factor(as.character(ddsHTSeq.mm.c5$treatment))
ddsHTSeq.mm.c5$treatment <- relevel(ddsHTSeq.mm.c5$treatment, ref = "WT.D0")
dds.c5 <- DESeq(ddsHTSeq.mm.c5)

res.c5 <- results(dds.c5, filterFun=ihw, lfcThreshold = 1)

summary(res.c5)

res.c5.plot <-  res.c5[!is.na(res.c5$padj),]

v.c5 <- data.frame(res.c5.plot[,c(2,6)])
colnames(v.c5)[2] <- "FDR"
ev.c5 <- EnhancedVolcano(v.c5,
                         lab = rownames(v.c5),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N56 D0 versus WT D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c5

res.sig.c5 <- res.c5[!is.na(res.c5$padj),]
res.sig.c5 <- res.sig.c5[abs(res.sig.c5$log2FoldChange) > 1,]
res.sig.c5 <- res.sig.c5[res.sig.c5$padj < 0.1,]
res.sig.c5 <- res.sig.c5[order(res.sig.c5$pvalue),]
res.c5.df <- data.frame(res.c5)
res.sig.c5.df <- data.frame(res.sig.c5)
WriteXLS(c("res.c5.df","res.sig.c5.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N56.D0_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c6 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D3","WT.D3")]
ddsHTSeq.mm.c6$treatment <- factor(as.character(ddsHTSeq.mm.c6$treatment))
ddsHTSeq.mm.c6$treatment <- relevel(ddsHTSeq.mm.c6$treatment, ref = "WT.D3")
dds.c6 <- DESeq(ddsHTSeq.mm.c6)

res.c6 <- results(dds.c6, filterFun=ihw, lfcThreshold = 1)
summary(res.c6)

res.c6.plot <-  res.c6[!is.na(res.c6$padj),]
v.c6 <- data.frame(res.c6.plot[,c(2,6)])
colnames(v.c6)[2] <- "FDR"
ev.c6 <- EnhancedVolcano(v.c6,
                         lab = rownames(v.c6),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N56 D3 versus WT D3',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))

ev.c6
res.sig.c6 <- res.c6[!is.na(res.c6$padj),]
res.sig.c6 <- res.sig.c6[abs(res.sig.c6$log2FoldChange) > 1,]
res.sig.c6 <- res.sig.c6[res.sig.c6$padj < 0.1,]
res.sig.c6 <- res.sig.c6[order(res.sig.c6$pvalue),]
res.c6.df <- data.frame(res.c6)
res.sig.c6.df <- data.frame(res.sig.c6)
WriteXLS(c("res.c6.df","res.sig.c6.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N56.D3_vs_WT.D3.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c7 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D7","WT.D7")]
ddsHTSeq.mm.c7$treatment <- factor(as.character(ddsHTSeq.mm.c7$treatment))
ddsHTSeq.mm.c7$treatment <- relevel(ddsHTSeq.mm.c7$treatment, ref = "WT.D7")
dds.c7 <- DESeq(ddsHTSeq.mm.c7)

res.c7 <- results(dds.c7, filterFun=ihw, lfcThreshold = 1)

summary(res.c7)

res.c7.plot <-  res.c7[!is.na(res.c7$padj),]

v.c7 <- data.frame(res.c7.plot[,c(2,6)])
colnames(v.c7)[2] <- "FDR"
ev.c7 <- EnhancedVolcano(v.c7,
                         lab = rownames(v.c7),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N56 D7 versus WT D7',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c7
res.sig.c7 <- res.c7[!is.na(res.c7$padj),]
res.sig.c7 <- res.sig.c7[abs(res.sig.c7$log2FoldChange) > 1,]
res.sig.c7 <- res.sig.c7[res.sig.c7$padj < 0.1,]
res.sig.c7 <- res.sig.c7[order(res.sig.c7$pvalue),]
res.c7.df <- data.frame(res.c7)
res.sig.c7.df <- data.frame(res.sig.c7)
WriteXLS(c("res.c7.df","res.sig.c7.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N56.D7_vs_WT.D7.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c8 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D10","WT.D10")]
ddsHTSeq.mm.c8$treatment <- factor(as.character(ddsHTSeq.mm.c8$treatment))
ddsHTSeq.mm.c8$treatment <- relevel(ddsHTSeq.mm.c8$treatment, ref = "WT.D10")
dds.c8 <- DESeq(ddsHTSeq.mm.c8)

res.c8 <- results(dds.c8, filterFun=ihw, lfcThreshold = 1)

summary(res.c8)

res.c8.plot <-  res.c8[!is.na(res.c8$padj),]

v.c8 <- data.frame(res.c8.plot[,c(2,6)])
colnames(v.c8)[2] <- "FDR"
ev.c8 <- EnhancedVolcano(v.c8,
                         lab = rownames(v.c8),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N56 D10 versus WT D10',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c8
es.sig.c8 <- res.c8[!is.na(res.c8$padj),]
res.sig.c8 <- res.sig.c8[abs(res.sig.c8$log2FoldChange) > 1,]
res.sig.c8 <- res.sig.c8[res.sig.c8$padj < 0.1,]
res.sig.c8 <- res.sig.c8[order(res.sig.c8$pvalue),]
res.c8.df <- data.frame(res.c8)
res.sig.c8.df <- data.frame(res.sig.c8)
WriteXLS(c("res.c8.df","res.sig.c8.df"),"Threshold_1_RNAseq_mESC_Epcam_KO_comp_1_N56.D10_vs_WT.D10.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

max.genes <- if(nrow(res.sig.c1) >= 50) 50 else nrow(res.sig.c1)
vst.c1 <- vst(dds.c1)
hmd <- assay(vst.c1[row.names(res.sig.c1)[1:max.genes],])
colnames(hmd) <- dds.c1$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D0 vs WT.D0")
max.genes <- if(nrow(res.sig.c1) >= 1000) 1000 else nrow(res.sig.c1)

vst.c1 <- vst(dds.c1)

hmd <- assay(vst.c1[row.names(res.sig.c1)[1:max.genes],])
colnames(hmd) <- dds.c1$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D0 vs WT.D0")
max.genes <- if(nrow(res.sig.c2) >= 50) 50 else nrow(res.sig.c2)

vst.c2 <- vst(dds.c2)

hmd <- assay(vst.c2[row.names(res.sig.c2)[1:max.genes],])
colnames(hmd) <- dds.c2$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D3 vs WT.D3")
max.genes <- if(nrow(res.sig.c3) >= 50) 50 else nrow(res.sig.c3)

vst.c3 <- vst(dds.c3)

hmd <- assay(vst.c3[row.names(res.sig.c3)[1:max.genes],])
colnames(hmd) <- dds.c3$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D7 vs WT.D7")
max.genes <- if(nrow(res.sig.c4) >= 50) 50 else nrow(res.sig.c4)

vst.c4 <- vst(dds.c4)

hmd <- assay(vst.c4[row.names(res.sig.c4)[1:max.genes],])
colnames(hmd) <- dds.c4$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D10 vs WT.D10")
max.genes <- if(nrow(res.sig.c5) >= 50) 50 else nrow(res.sig.c5)

vst.c5 <- vst(dds.c5)

hmd <- assay(vst.c5[row.names(res.sig.c5)[1:max.genes],])
colnames(hmd) <- dds.c5$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D0 vs WT.D0")
max.genes <- if(nrow(res.sig.c6) >= 50) 50 else nrow(res.sig.c6)

vst.c6 <- vst(dds.c6)

hmd <- assay(vst.c6[row.names(res.sig.c6)[1:max.genes],])
colnames(hmd) <- dds.c6$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D3 vs WT.D3")
max.genes <- if(nrow(res.sig.c7) >= 50) 50 else nrow(res.sig.c7)

vst.c7 <- vst(dds.c7)

hmd <- assay(vst.c7[row.names(res.sig.c7)[1:max.genes],])
colnames(hmd) <- dds.c7$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D7 vs WT.D7")
max.genes <- if(nrow(res.sig.c4) >= 50) 50 else nrow(res.sig.c4)

st.c4 <- vst(dds.c4)

hmd <- assay(vst.c4[row.names(res.sig.c4)[1:max.genes],])
colnames(hmd) <- dds.c4$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D10 vs WT.D10")
max.genes <- if(nrow(res.sig.c6) >= 50) 50 else nrow(res.sig.c6)

vst.c6 <- vst(dds.c6)

hmd <- assay(vst.c6[row.names(res.sig.c6)[1:max.genes],])
colnames(hmd) <- dds.c6$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D3 vs WT.D3")
max.genes <- if(nrow(res.sig.c7) >= 50) 50 else nrow(res.sig.c7)

vst.c7 <- vst(dds.c7)

hmd <- assay(vst.c7[row.names(res.sig.c7)[1:max.genes],])
colnames(hmd) <- dds.c7$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D7 vs WT.D7")

max.genes <- if(nrow(res.sig.c8) >= 50) 50 else nrow(res.sig.c8)

vst.c8 <- vst(dds.c8)

hmd <- assay(vst.c8[row.names(res.sig.c8)[1:max.genes],])
colnames(hmd) <- dds.c8$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D10 vs WT.D10")


## ---- DE along time (design ~ time + treat; c3/c4 re-defined) -------------
##overall
ddsHTSeq.mm$treatment

time.vec <- ddsHTSeq.mm$treatment
time.vec <- gsub("WT.","",time.vec)
time.vec <- gsub("N56.","",time.vec)
time.vec <- gsub("N114.","",time.vec)
time.vec <- factor(time.vec)

treat.vec <- ddsHTSeq.mm$treatment
treat.vec <- gsub(".D0","",treat.vec)
treat.vec <- gsub(".D3","",treat.vec)
treat.vec <- gsub(".D7","",treat.vec)
treat.vec <- gsub(".D10","",treat.vec)
treat.vec <- factor(treat.vec)
treat.vec <- relevel(treat.vec, ref = "WT")

ddsHTSeq.mm$time <- time.vec
ddsHTSeq.mm$treat <- treat.vec

ddsHTSeq.mm.c3 <- ddsHTSeq.mm[,ddsHTSeq.mm$treat%in%c("WT","N114")]

ddsHTSeq.mm.c3 <- DESeqDataSetFromMatrix(counts(ddsHTSeq.mm.c3), colData = colData(ddsHTSeq.mm.c3), design = time~treat)

dds.c3 <- DESeq(ddsHTSeq.mm.c3)

res.c3 <- results(dds.c3, filterFun=ihw, lfcThreshold = 0.5)

summary(res.c3)

res.c3.plot <-  res.c3[!is.na(res.c3$padj),]


res.sig.c3 <- res.c3[!is.na(res.c3$padj),]
res.sig.c3 <- res.sig.c3[abs(res.sig.c3$log2FoldChange) > 1,]
res.sig.c3 <- res.sig.c3[res.sig.c3$padj < 0.1,]
res.sig.c3 <- res.sig.c3[order(res.sig.c3$pvalue),]

max.genes <- if(nrow(res.sig.c3) >= 50) 50 else nrow(res.sig.c3)

vst.c3 <- vst(dds.c3)

hmd <- assay(vst.c3[row.names(res.sig.c3)[1:max.genes],])
colnames(hmd) <- dds.c3$treatment

hmd <- hmd[,mixedsort(colnames(hmd))]

pheatmap(hmd, col = bluered(50), cluster_rows = T, cluster_cols = F, scale = "row", main = "Dex genes N56 vs WT")

ddsHTSeq.mm.c4 <- ddsHTSeq.mm[,ddsHTSeq.mm$treat%in%c("WT","N56")]

ddsHTSeq.mm.c4 <- DESeqDataSetFromMatrix(counts(ddsHTSeq.mm.c4), colData = colData(ddsHTSeq.mm.c4), design = time~treat)

dds.c4 <- DESeq(ddsHTSeq.mm.c4)

res.c4 <- results(dds.c4, filterFun=ihw, lfcThreshold = 1)

summary(res.c4)

res.c4.plot <-  res.c4[!is.na(res.c4$padj),]


res.sig.c4 <- res.c4[!is.na(res.c4$padj),]
res.sig.c4 <- res.sig.c4[abs(res.sig.c4$log2FoldChange) > 1,]
res.sig.c4 <- res.sig.c4[res.sig.c4$padj < 0.1,]
res.sig.c4 <- res.sig.c4[order(res.sig.c4$pvalue),]

max.genes <- if(nrow(res.sig.c4) >= 50) 50 else nrow(res.sig.c4)

vst.c4 <- vst(dds.c4)

hmd <- assay(vst.c4[row.names(res.sig.c4)[1:max.genes],])
colnames(hmd) <- dds.c4$treatment
hmd <- hmd[,mixedsort(colnames(hmd))]

pheatmap(hmd, col = bluered(50), cluster_rows = T, cluster_cols = F, scale = "row", main = "Dex genes N56 vs WT")


## ---- DE along differentiation within each genotype vs D0 (c9-c17) --------
ddsHTSeq.mm.c9 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("WT.D3","WT.D0")]
ddsHTSeq.mm.c9$treatment <- factor(as.character(ddsHTSeq.mm.c9$treatment))
ddsHTSeq.mm.c9$treatment <- relevel(ddsHTSeq.mm.c9$treatment, ref = "WT.D0")
dds.c9 <- DESeq(ddsHTSeq.mm.c9)

res.c9 <- results(dds.c9, filterFun=ihw, lfcThreshold = 1)

summary(res.c9)

res.c9.plot <-  res.c9[!is.na(res.c9$padj),]

v.c9 <- data.frame(res.c9.plot[,c(2,6)])
colnames(v.c9)[2] <- "FDR"
ev.c9 <- EnhancedVolcano(v.c9,
                         lab = rownames(v.c9),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'WT D3 versus WT D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c9
es.sig.c9 <- res.c9[!is.na(res.c9$padj),]
res.sig.c9 <- res.sig.c9[abs(res.sig.c9$log2FoldChange) > 1,]
res.sig.c9 <- res.sig.c9[res.sig.c9$padj < 0.1,]
res.sig.c9 <- res.sig.c9[order(res.sig.c9$pvalue),]
res.c9.df <- data.frame(res.c9)
res.sig.c9.df <- data.frame(res.sig.c9)
WriteXLS(c("res.c9.df","res.sig.c9.df"),"Threshold_1_RNAseq_mESC_Epcam_WT.D3_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


ddsHTSeq.mm.c10 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("WT.D7","WT.D0")]
ddsHTSeq.mm.c10$treatment <- factor(as.character(ddsHTSeq.mm.c10$treatment))
ddsHTSeq.mm.c10$treatment <- relevel(ddsHTSeq.mm.c9$treatment, ref = "WT.D0")
dds.c10 <- DESeq(ddsHTSeq.mm.c10)

res.c10 <- results(dds.c10, filterFun=ihw, lfcThreshold = 1)

summary(res.c10)

res.c10.plot <-  res.c10[!is.na(res.c10$padj),]

v.c10 <- data.frame(res.c10.plot[,c(2,6)])
colnames(v.c10)[2] <- "FDR"
ev.c10 <- EnhancedVolcano(v.c10,
                         lab = rownames(v.c10),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'WT D7 versus WT D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c10
res.sig.c10 <- res.c10[!is.na(res.c10$padj),]
res.sig.c10 <- res.sig.c10[abs(res.sig.c10$log2FoldChange) > 1,]
res.sig.c10 <- res.sig.c10[res.sig.c10$padj < 0.1,]
res.sig.c10 <- res.sig.c10[order(res.sig.c10$pvalue),]
res.c10.df <- data.frame(res.c10)
res.sig.c10.df <- data.frame(res.sig.c10)
WriteXLS(c("res.c10.df","res.sig.c10.df"),"Threshold_1_RNAseq_mESC_Epcam_WT.D7_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c11 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("WT.D10","WT.D0")]
ddsHTSeq.mm.c11$treatment <- factor(as.character(ddsHTSeq.mm.c11$treatment))
ddsHTSeq.mm.c11$treatment <- relevel(ddsHTSeq.mm.c11$treatment, ref = "WT.D0")
dds.c11 <- DESeq(ddsHTSeq.mm.c11)

res.c11 <- results(dds.c11, filterFun=ihw, lfcThreshold = 1)

summary(res.c11)

res.c11.plot <-  res.c11[!is.na(res.c11$padj),]

v.c11 <- data.frame(res.c11.plot[,c(2,6)])
colnames(v.c11)[2] <- "FDR"
ev.c11 <- EnhancedVolcano(v.c11,
                         lab = rownames(v.c11),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'WT D10 versus WT D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c11
res.sig.c11 <- res.c11[!is.na(res.c11$padj),]
res.sig.c11 <- res.sig.c11[abs(res.sig.c11$log2FoldChange) > 1,]
res.sig.c11 <- res.sig.c11[res.sig.c11$padj < 0.1,]
res.sig.c11 <- res.sig.c11[order(res.sig.c11$pvalue),]
res.c11.df <- data.frame(res.c11)
res.sig.c11.df <- data.frame(res.sig.c11)
WriteXLS(c("res.c11.df","res.sig.c11.df"),"Threshold_1_RNAseq_mESC_Epcam_WT.D10_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


ddsHTSeq.mm.c12 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","N114.D0")]
ddsHTSeq.mm.c12$treatment <- factor(as.character(ddsHTSeq.mm.c12$treatment))
ddsHTSeq.mm.c12$treatment <- relevel(ddsHTSeq.mm.c12$treatment, ref = "N114.D0")
dds.c12 <- DESeq(ddsHTSeq.mm.c12)

res.c12 <- results(dds.c12, filterFun=ihw, lfcThreshold = 1)

summary(res.c12)

res.c12.plot <-  res.c12[!is.na(res.c12$padj),]

v.c12 <- data.frame(res.c12.plot[,c(2,6)])
colnames(v.c12)[2] <- "FDR"
ev.c12 <- EnhancedVolcano(v.c12,
                         lab = rownames(v.c12),
                         x = 'log2FoldChange',
                         y = 'FDR',
                         title = 'N114 D3 versus N114 D0',
                         pCutoff = 0.25,
                         FCcutoff = 0,
                         pointSize = 3.0,
                         labSize = 3.0,
                         ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c12
res.sig.c12 <- res.c12[!is.na(res.c12$padj),]
res.sig.c12 <- res.sig.c12[abs(res.sig.c12$log2FoldChange) > 1,]
res.sig.c12 <- res.sig.c12[res.sig.c12$padj < 0.1,]
res.sig.c12 <- res.sig.c12[order(res.sig.c12$pvalue),]
res.c12.df <- data.frame(res.c12)
res.sig.c12.df <- data.frame(res.sig.c12)
WriteXLS(c("res.c12.df","res.sig.c12.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D3_vs_N114.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c13 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D7","N114.D0")]
ddsHTSeq.mm.c13$treatment <- factor(as.character(ddsHTSeq.mm.c13$treatment))
ddsHTSeq.mm.c13$treatment <- relevel(ddsHTSeq.mm.c13$treatment, ref = "N114.D0")
dds.c13 <- DESeq(ddsHTSeq.mm.c13)

res.c13 <- results(dds.c13, filterFun=ihw, lfcThreshold = 1)

summary(res.c13)

res.c13.plot <-  res.c13[!is.na(res.c13$padj),]

v.c13 <- data.frame(res.c13.plot[,c(2,6)])
colnames(v.c13)[2] <- "FDR"
ev.c13 <- EnhancedVolcano(v.c13,
                          lab = rownames(v.c13),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D7 versus N114 D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c13
res.sig.c13 <- res.c13[!is.na(res.c13$padj),]
res.sig.c13 <- res.sig.c13[abs(res.sig.c13$log2FoldChange) > 1,]
res.sig.c13 <- res.sig.c13[res.sig.c13$padj < 0.1,]
res.sig.c13 <- res.sig.c13[order(res.sig.c13$pvalue),]
res.c13.df <- data.frame(res.c13)
res.sig.c13.df <- data.frame(res.sig.c13)
WriteXLS(c("res.c13.df","res.sig.c13.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D7_vs_N114.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)

ddsHTSeq.mm.c14 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D10","N114.D0")]
ddsHTSeq.mm.c14$treatment <- factor(as.character(ddsHTSeq.mm.c14$treatment))
ddsHTSeq.mm.c14$treatment <- relevel(ddsHTSeq.mm.c14$treatment, ref = "N114.D0")
dds.c14 <- DESeq(ddsHTSeq.mm.c14)

res.c14 <- results(dds.c14, filterFun=ihw, lfcThreshold = 1)

summary(res.c14)

res.c14.plot <-  res.c14[!is.na(res.c14$padj),]

v.c14 <- data.frame(res.c14.plot[,c(2,6)])
colnames(v.c14)[2] <- "FDR"
ev.c14 <- EnhancedVolcano(v.c14,
                          lab = rownames(v.c14),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D10 versus N114 D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c14
res.sig.c14 <- res.c14[!is.na(res.c14$padj),]
res.sig.c14 <- res.sig.c14[abs(res.sig.c14$log2FoldChange) > 1,]
res.sig.c14 <- res.sig.c14[res.sig.c14$padj < 0.1,]
res.sig.c14 <- res.sig.c14[order(res.sig.c14$pvalue),]
res.c14.df <- data.frame(res.c14)
res.sig.c14.df <- data.frame(res.sig.c14)
WriteXLS(c("res.c14.df","res.sig.c14.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D10_vs_N114.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


ddsHTSeq.mm.c15 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D3","N56.D0")]
ddsHTSeq.mm.c15$treatment <- factor(as.character(ddsHTSeq.mm.c15$treatment))
ddsHTSeq.mm.c15$treatment <- relevel(ddsHTSeq.mm.c15$treatment, ref = "N56.D0")
dds.c15 <- DESeq(ddsHTSeq.mm.c15)

res.c15 <- results(dds.c15, filterFun=ihw, lfcThreshold = 1)

summary(res.c15)

res.c15.plot <-  res.c15[!is.na(res.c15$padj),]

v.c15 <- data.frame(res.c15.plot[,c(2,6)])
colnames(v.c15)[2] <- "FDR"
ev.c15 <- EnhancedVolcano(v.c15,
                          lab = rownames(v.c15),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N56 D3 versus N56 D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c15
res.sig.c15 <- res.c15[!is.na(res.c15$padj),]
res.sig.c15 <- res.sig.c15[abs(res.sig.c15$log2FoldChange) > 1,]
res.sig.c15 <- res.sig.c15[res.sig.c15$padj < 0.1,]
res.sig.c15 <- res.sig.c15[order(res.sig.c15$pvalue),]
res.c15.df <- data.frame(res.c15)
res.sig.c15.df <- data.frame(res.sig.c15)
WriteXLS(c("res.c15.df","res.sig.c15.df"),"Threshold_1_RNAseq_mESC_Epcam_N56.D3_vs_N56.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


ddsHTSeq.mm.c16 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D7","N56.D0")]
ddsHTSeq.mm.c16$treatment <- factor(as.character(ddsHTSeq.mm.c16$treatment))
ddsHTSeq.mm.c16$treatment <- relevel(ddsHTSeq.mm.c16$treatment, ref = "N56.D0")
dds.c16 <- DESeq(ddsHTSeq.mm.c16)

res.c16 <- results(dds.c16, filterFun=ihw, lfcThreshold = 1)

summary(res.c16)

res.c16.plot <-  res.c16[!is.na(res.c16$padj),]

v.c16 <- data.frame(res.c16.plot[,c(2,6)])
colnames(v.c16)[2] <- "FDR"
ev.c16 <- EnhancedVolcano(v.c16,
                          lab = rownames(v.c16),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N56 D7 versus N56 D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c16
res.sig.c16 <- res.c16[!is.na(res.c16$padj),]
res.sig.c16 <- res.sig.c16[abs(res.sig.c16$log2FoldChange) > 1,]
res.sig.c16 <- res.sig.c16[res.sig.c16$padj < 0.1,]
res.sig.c16 <- res.sig.c16[order(res.sig.c16$pvalue),]
res.c16.df <- data.frame(res.c16)
res.sig.c16.df <- data.frame(res.sig.c16)
WriteXLS(c("res.c16.df","res.sig.c16.df"),"Threshold_1_RNAseq_mESC_Epcam_N56.D7_vs_N56.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


ddsHTSeq.mm.c17 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N56.D10","N56.D0")]
ddsHTSeq.mm.c17$treatment <- factor(as.character(ddsHTSeq.mm.c17$treatment))
ddsHTSeq.mm.c17$treatment <- relevel(ddsHTSeq.mm.c17$treatment, ref = "N56.D0")
dds.c17 <- DESeq(ddsHTSeq.mm.c17)

res.c17 <- results(dds.c17, filterFun=ihw, lfcThreshold = 1)

summary(res.c17)

res.c17.plot <-  res.c17[!is.na(res.c17$padj),]

v.c17 <- data.frame(res.c17.plot[,c(2,6)])
colnames(v.c17)[2] <- "FDR"
ev.c17 <- EnhancedVolcano(v.c17,
                          lab = rownames(v.c17),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N56 D10 versus N56 D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c17
res.sig.c17 <- res.c17[!is.na(res.c17$padj),]
res.sig.c17 <- res.sig.c17[abs(res.sig.c17$log2FoldChange) > 1,]
res.sig.c17 <- res.sig.c17[res.sig.c17$padj < 0.1,]
res.sig.c17 <- res.sig.c17[order(res.sig.c17$pvalue),]
res.c17.df <- data.frame(res.c17)
res.sig.c17.df <- data.frame(res.sig.c17)
WriteXLS(c("res.c17.df","res.sig.c17.df"),"Threshold_1_RNAseq_mESC_Epcam_N56.D10_vs_N56.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)


## ---- Pooled KO (#56 + #114) vs WT per time point (c18-c21) --------------
ddsHTSeq.mm.c18 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D0","N56.D0", "WT.D0")]
ddsHTSeq.mm.c18$treatment <- factor(as.character(ddsHTSeq.mm.c18$treatment))
ddsHTSeq.mm.c18$treatment <- relevel(ddsHTSeq.mm.c18$treatment, ref = "WT.D0")
dds.c18 <- DESeq(ddsHTSeq.mm.c18)

res.c18 <- results(dds.c18, filterFun=ihw, lfcThreshold = 1)

summary(res.c18)

res.c18.plot <-  res.c18[!is.na(res.c18$padj),]

v.c18 <- data.frame(res.c18.plot[,c(2,6)])
colnames(v.c18)[2] <- "FDR"
# Volcano plot (larger points/labels)
ev.c18 <- EnhancedVolcano(v.c18,
                          lab = rownames(v.c18),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D0 N56 D0 versus WT D0',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 6.0,
                          labSize = 6.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c18


res.sig.c18 <- res.c18[!is.na(res.c18$padj),]
res.sig.c18 <- res.sig.c18[abs(res.sig.c18$log2FoldChange) > 1,]
res.sig.c18 <- res.sig.c18[res.sig.c18$padj < 0.1,]
res.sig.c18 <- res.sig.c18[order(res.sig.c18$pvalue),]
res.c18.df <- data.frame(res.c18)
res.sig.c18.df <- data.frame(res.sig.c18)
WriteXLS(c("res.c18.df","res.sig.c18.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D0_N56.D0_vs_WT.D0.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)
# Find up-regulated genes
up_regulated_genes_D0 <- res.sig.c18.df[res.sig.c18.df$log2FoldChange > 0, ]

# Find down-regulated genes
down_regulated_genes_D0 <- res.sig.c18.df[res.sig.c18.df$log2FoldChange < 0, ]

n_up_regulated_genes_D0 <- nrow(up_regulated_genes_D0)
n_down_regulated_genes_D0 <- nrow(down_regulated_genes_D0)
n_up_regulated_genes_D0
n_down_regulated_genes_D0


ddsHTSeq.mm.c19 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","N56.D3", "WT.D3")]
ddsHTSeq.mm.c19$treatment <- factor(as.character(ddsHTSeq.mm.c19$treatment))
ddsHTSeq.mm.c19$treatment <- relevel(ddsHTSeq.mm.c19$treatment, ref = "WT.D3")
dds.c19 <- DESeq(ddsHTSeq.mm.c19)

res.c19 <- results(dds.c19, filterFun=ihw, lfcThreshold = 1)

summary(res.c19)

res.c19.plot <-  res.c19[!is.na(res.c19$padj),]

v.c19 <- data.frame(res.c19.plot[,c(2,6)])
colnames(v.c19)[2] <- "FDR"
# Volcano plot (larger points/labels)
ev.c19 <- EnhancedVolcano(v.c19,
                          lab = rownames(v.c19),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D3 N56 D3 versus WT D3',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 6.0,
                          labSize = 6.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c19


res.sig.c19 <- res.c19[!is.na(res.c19$padj),]
res.sig.c19 <- res.sig.c19[abs(res.sig.c19$log2FoldChange) > 1,]
res.sig.c19 <- res.sig.c19[res.sig.c19$padj < 0.1,]
res.sig.c19 <- res.sig.c19[order(res.sig.c19$pvalue),]
res.c19.df <- data.frame(res.c19)
res.sig.c19.df <- data.frame(res.sig.c19)
WriteXLS(c("res.c19.df","res.sig.c19.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D3_N56.D3_vs_WT.D3.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)
# Find up-regulated genes
up_regulated_genes_D3 <- res.sig.c19.df[res.sig.c19.df$log2FoldChange > 0, ]

# Find down-regulated genes
down_regulated_genes_D3 <- res.sig.c19.df[res.sig.c19.df$log2FoldChange < 0, ]

n_up_regulated_genes_D3 <- nrow(up_regulated_genes_D3)
n_down_regulated_genes_D3 <- nrow(down_regulated_genes_D3)
n_up_regulated_genes_D3
n_down_regulated_genes_D3


ddsHTSeq.mm.c20 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D7","N56.D7", "WT.D7")]
ddsHTSeq.mm.c20$treatment <- factor(as.character(ddsHTSeq.mm.c20$treatment))
ddsHTSeq.mm.c20$treatment <- relevel(ddsHTSeq.mm.c20$treatment, ref = "WT.D7")
dds.c20 <- DESeq(ddsHTSeq.mm.c20)

res.c20 <- results(dds.c20, filterFun=ihw, lfcThreshold = 1)

summary(res.c20)

res.c20.plot <-  res.c20[!is.na(res.c20$padj),]

v.c20 <- data.frame(res.c20.plot[,c(2,6)])
colnames(v.c20)[2] <- "FDR"
# Volcano plot (larger points/labels)
ev.c20 <- EnhancedVolcano(v.c20,
                          lab = rownames(v.c20),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D7 N56 D7 versus WT D7',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 6.0,
                          labSize = 6.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c20


res.sig.c20 <- res.c20[!is.na(res.c20$padj),]
res.sig.c20 <- res.sig.c20[abs(res.sig.c20$log2FoldChange) > 1,]
res.sig.c20 <- res.sig.c20[res.sig.c20$padj < 0.1,]
res.sig.c20 <- res.sig.c20[order(res.sig.c20$pvalue),]
res.c20.df <- data.frame(res.c20)
res.sig.c20.df <- data.frame(res.sig.c20)
WriteXLS(c("res.c20.df","res.sig.c20.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D7_N56.D7_vs_WT.D7.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)
# Find up-regulated genes
up_regulated_genes_D7 <- res.sig.c20.df[res.sig.c20.df$log2FoldChange > 0, ]

# Find down-regulated genes
down_regulated_genes_D7 <- res.sig.c20.df[res.sig.c20.df$log2FoldChange < 0, ]

n_up_regulated_genes_D7 <- nrow(up_regulated_genes_D7)
n_down_regulated_genes_D7 <- nrow(down_regulated_genes_D7)
n_up_regulated_genes_D7
n_down_regulated_genes_D7


ddsHTSeq.mm.c21 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D10","N56.D10", "WT.D10")]
ddsHTSeq.mm.c21$treatment <- factor(as.character(ddsHTSeq.mm.c21$treatment))
ddsHTSeq.mm.c21$treatment <- relevel(ddsHTSeq.mm.c21$treatment, ref = "WT.D10")
dds.c21 <- DESeq(ddsHTSeq.mm.c21)

res.c21 <- results(dds.c21, filterFun=ihw, lfcThreshold = 1)

summary(res.c21)

res.c21.plot <-  res.c21[!is.na(res.c21$padj),]

v.c21 <- data.frame(res.c21.plot[,c(2,6)])
colnames(v.c21)[2] <- "FDR"
# Volcano plot (larger points/labels)
ev.c21 <- EnhancedVolcano(v.c21,
                          lab = rownames(v.c21),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D10 N56 D10 versus WT D10',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 6.0,
                          labSize = 6.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c21


res.sig.c21 <- res.c21[!is.na(res.c21$padj),]
res.sig.c21 <- res.sig.c21[abs(res.sig.c21$log2FoldChange) > 1,]
res.sig.c21 <- res.sig.c21[res.sig.c21$padj < 0.1,]
res.sig.c21 <- res.sig.c21[order(res.sig.c21$pvalue),]
res.c21.df <- data.frame(res.c21)
res.sig.c21.df <- data.frame(res.sig.c21)
WriteXLS(c("res.c21.df","res.sig.c21.df"),"Threshold_1_RNAseq_mESC_Epcam_N114.D10_N56.D10_vs_WT.D10.xlsx", SheetNames = c("all.results","stat.sig.results"), row.names = T)
# Find up-regulated genes
up_regulated_genes_D10 <- res.sig.c21.df[res.sig.c21.df$log2FoldChange > 0, ]

# Find down-regulated genes
down_regulated_genes_D10 <- res.sig.c21.df[res.sig.c21.df$log2FoldChange < 0, ]

n_up_regulated_genes_D10 <- nrow(up_regulated_genes_D10)
n_down_regulated_genes_D10 <- nrow(down_regulated_genes_D10)
n_up_regulated_genes_D10
n_down_regulated_genes_D10

wb <- createWorkbook()
addWorksheet(wb, "up_regulated_genes")
writeData(wb, "up_regulated_genes", up_regulated_genes_D10)
saveWorkbook(wb, "D10_up_regulated_genes.xlsx", overwrite = TRUE)


## ---- Heatmaps of top-50 DEGs for each contrast ---------------------------
max.genes <- if(nrow(res.sig.c9) >= 50) 50 else nrow(res.sig.c9)

vst.c9 <- vst(dds.c9)

hmd <- assay(vst.c9[row.names(res.sig.c9)[1:max.genes],])
colnames(hmd) <- dds.c9$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes WT.D3 vs WT.D0")


max.genes <- if(nrow(res.sig.c10) >= 50) 50 else nrow(res.sig.c10)

vst.c10 <- vst(dds.c10)

hmd <- assay(vst.c10[row.names(res.sig.c10)[1:max.genes],])
colnames(hmd) <- dds.c10$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes WT.D7 vs WT.D0")


max.genes <- if(nrow(res.sig.c11) >= 50) 50 else nrow(res.sig.c11)

vst.c11 <- vst(dds.c11)

hmd <- assay(vst.c11[row.names(res.sig.c11)[1:max.genes],])
colnames(hmd) <- dds.c11$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes WT.D10 vs WT.D0")


max.genes <- if(nrow(res.sig.c12) >= 50) 50 else nrow(res.sig.c12)

vst.c12 <- vst(dds.c12)

hmd <- assay(vst.c12[row.names(res.sig.c12)[1:max.genes],])
colnames(hmd) <- dds.c12$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D3 vs N114.D0")


max.genes <- if(nrow(res.sig.c13) >= 50) 50 else nrow(res.sig.c13)

vst.c13 <- vst(dds.c13)

hmd <- assay(vst.c13[row.names(res.sig.c13)[1:max.genes],])
colnames(hmd) <- dds.c13$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D7 vs N114.D0")


max.genes <- if(nrow(res.sig.c14) >= 50) 50 else nrow(res.sig.c14)

vst.c14 <- vst(dds.c14)

hmd <- assay(vst.c14[row.names(res.sig.c14)[1:max.genes],])
colnames(hmd) <- dds.c14$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D10 vs N114.D0")


max.genes <- if(nrow(res.sig.c15) >= 50) 50 else nrow(res.sig.c15)

vst.c15 <- vst(dds.c15)

hmd <- assay(vst.c15[row.names(res.sig.c15)[1:max.genes],])
colnames(hmd) <- dds.c15$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D3 vs N56.D0")


max.genes <- if(nrow(res.sig.c16) >= 50) 50 else nrow(res.sig.c16)

vst.c16 <- vst(dds.c16)

hmd <- assay(vst.c16[row.names(res.sig.c16)[1:max.genes],])
colnames(hmd) <- dds.c16$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D7 vs N56.D0")


max.genes <- if(nrow(res.sig.c17) >= 50) 50 else nrow(res.sig.c17)

vst.c17 <- vst(dds.c17)

hmd <- assay(vst.c17[row.names(res.sig.c17)[1:max.genes],])
colnames(hmd) <- dds.c17$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N56.D10 vs N56.D0")


max.genes <- if(nrow(res.sig.c18) >= 50) 50 else nrow(res.sig.c18)

vst.c18 <- vst(dds.c18)

hmd <- assay(vst.c18[row.names(res.sig.c18)[1:max.genes],])
colnames(hmd) <- dds.c18$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D0 N56.D0 vs WT.D0")


max.genes <- if(nrow(res.sig.c19) >= 50) 50 else nrow(res.sig.c19)

vst.c19 <- vst(dds.c19)

hmd <- assay(vst.c19[row.names(res.sig.c19)[1:max.genes],])
colnames(hmd) <- dds.c19$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D3 N56.D3 vs WT.D3")


max.genes <- if(nrow(res.sig.c20) >= 50) 50 else nrow(res.sig.c20)

vst.c20 <- vst(dds.c20)

hmd <- assay(vst.c20[row.names(res.sig.c20)[1:max.genes],])
colnames(hmd) <- dds.c20$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D7 N56.D7 vs WT.D7")


max.genes <- if(nrow(res.sig.c21) >= 50) 50 else nrow(res.sig.c21)

vst.c21 <- vst(dds.c21)

hmd <- assay(vst.c21[row.names(res.sig.c21)[1:max.genes],])
colnames(hmd) <- dds.c21$treatment

pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Dex genes N114.D10 N56.D10 vs WT.D10")


## ---- GO over-representation (D3 pooled KO vs WT) ------------------------
rld.c19 <- rlog(dds.c19)
res.c19 <- results(dds.c19, filterFun=ihw, lfcThreshold = 1)
results <- results(dds.c19)
significant_genes <- as.data.frame(results) %>%  # convert results object to a data frame
  mutate(gene = row.names(.)) %>%  # create a new column "gene" with the row names (i.e., gene names)
  filter(padj < 0.05) %>%
  select(gene, log2FoldChange, padj)
ora_results <- enrichGO(significant_genes$gene,
                        universe = names(org.Mm.eg.db),
                        OrgDb = org.Mm.eg.db,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 5,
                        maxGSSize = 500)


## ---- Target-gene heatmaps / violin plots (D0, D3, D7, D10 Venn genes) ----
## D0 Venn genes heat maps

genes_of_interest <- c("Nr0b1", "Inhbb", "Gata6", "Foxc1", "Nanog", "Esrrb", "Stat3", "Foxa2", "Kdm6a", "HMGa2", "Col5a2", "Col5a1")
vst.subset <- vst.c18[row.names(vst.c18) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c18$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D0 Target genes")

## CD86 target gene expression Violin plot D0
# Define your data
genes_of_interest <- c("Cd86")
vst.subset <- vst.c18[row.names(vst.c18) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c18$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Cd86" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Cd86" = "black", "other genes" = "black")) +
  labs(title = "CD86 expression D0", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))


## CD86 target gene expression Violin plot D3
genes_of_interest <- c("Cd86")
vst.subset <- vst.c19[row.names(vst.c19) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Cd86" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Cd86" = "black", "other genes" = "black")) +
  labs(title = "CD86 expression D3", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

## CD86 target gene expression Violin plot D7
genes_of_interest <- c("Cd86")
vst.subset <- vst.c20[row.names(vst.c20) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c20$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Cd86" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Cd86" = "black", "other genes" = "black")) +
  labs(title = "CD86 expression D7", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

## CD86 target gene expression Violin plot D10
genes_of_interest <- c("Cd86")
vst.subset <- vst.c21[row.names(vst.c21) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c21$treatment

## Hand1
## Hand1 target gene expression Violin plot D0
# Define your data
genes_of_interest <- c("Hand1")
vst.subset <- vst.c18[row.names(vst.c18) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c18$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Hand1" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Hand1" = "black", "other genes" = "black")) +
  labs(title = "Hand1 expression D0", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))


## Hand1 target gene expression Violin plot D3
genes_of_interest <- c("Hand1")
vst.subset <- vst.c19[row.names(vst.c19) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Hand1" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Hand1" = "black", "other genes" = "black")) +
  labs(title = "Hand1 expression D3", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

## Hand1 target gene expression Violin plot D7
genes_of_interest <- c("Hand1")
vst.subset <- vst.c20[row.names(vst.c20) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c20$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Hand1" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Hand1" = "black", "other genes" = "black")) +
  labs(title = "Hand1 expression D7", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

## Hand1 target gene expression Violin plot D10
genes_of_interest <- c("Hand1")
vst.subset <- vst.c21[row.names(vst.c21) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c21$treatment

# Convert the data to long format
df <- data.frame(gene = rownames(hmd),
                 sample = colnames(hmd),
                 expression = as.vector(hmd))

# Create the violin plot with replicates using ggplot2
ggplot(df, aes(x = sample, y = expression, fill = gene, group = sample)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_point(aes(color = gene), size = 2, position = position_jitter(width = 0.2)) +
  scale_fill_manual(values = c("Hand1" = "magenta", "other genes" = "black")) +
  scale_color_manual(values = c("Hand1" = "black", "other genes" = "black")) +
  labs(title = "Hand1 expression D10", x = "Sample", y = "Expression") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22))

for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Print the gene name and mean values
  print(paste(gene,":", sep = " "))
  print(mean_values)
}

mean_values_list <- list()
for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Add the gene name and mean values to a list
  mean_values_list <- rbind(mean_values_list, c(gene,mean_values))
}
# Write the list to a csv file
write.csv(mean_values_list, file = "mean_values.csv", row.names = FALSE)


## D3 Venn genes heat maps

genes_of_interest <- c("Gata6", "Cdkn2a", "T", "Dkk1", "Sall3", "Foxa2", "Wnt11", "Mesp1", "Eomes", "Mixl1", "Fgf8", "Fgf10", "Uty")
vst.subset <- vst.c19[row.names(vst.c19) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D3 Target genes")

print(hmd)
wb <- createWorkbook()
addWorksheet(wb, "hmd")
writeData(wb, "hmd", hmd)
saveWorkbook(wb, "hmd.xlsx", overwrite = TRUE)


for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Print the gene name and mean values
  print(paste(gene,":", sep = " "))
  print(mean_values)
}

mean_values_list <- list()
for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Add the gene name and mean values to a list
  mean_values_list <- rbind(mean_values_list, c(gene,mean_values))
}
# Write the list to a csv file
write.csv(mean_values_list, file = "D3_mean_values.csv", row.names = FALSE)


## D7 Venn genes heat maps
genes_of_interest <- c("Esrrb", "Nanog", "Utf1", "Zic3", "Fgf4", "Lefty1", "Twist1", "Fbxo15", "Nr5a2", "Klf5", "Tfcp2l1", "Foxc1", "Foxa2", "Col5a2", "Col5a1")
vst.subset <- vst.c20[row.names(vst.c20) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c20$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D7 Target genes")

print(hmd)
wb <- createWorkbook()
addWorksheet(wb, "hmd")
writeData(wb, "hmd", hmd)
saveWorkbook(wb, "hmd.xlsx", overwrite = TRUE)


mean_values_list <- list()
for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Add the gene name and mean values to a list
  mean_values_list <- rbind(mean_values_list, c(gene,mean_values))
}
# Write the list to a csv file
write.csv(mean_values_list, file = "D7_mean_values.csv", row.names = FALSE)


## D10 Venn genes

genes_of_interest <- c("Esrrb", "Nanog","Dazl", "Fbxo15", "Cdx2", "Hand1", "Fgf4")
vst.subset <- vst.c21[row.names(vst.c21) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c21$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D10 Target genes")

print(hmd)
wb <- createWorkbook()
addWorksheet(wb, "hmd")
writeData(wb, "hmd", hmd)
saveWorkbook(wb, "hmd.xlsx", overwrite = TRUE)


## ---- D3 target genes: replicate-level tables and summaries --------------
colData(dds.c19)$replicate <- c("WT_D3_1", "WT_D3_2", "WT_D3_3", "WT_D3_4", "N114_D3_1", "N114_D3_2", "N114_D3_3", "N114_D3_4", "N56_D3_1", "N56_D3_2", "N56_D3_3", "N56_D3_4")
genes_of_interest <- c("Gata6", "Cdkn2a", "T", "Dkk1", "Sall3", "Foxa2", "Wnt11", "Mesp1", "Eomes", "Mixl1", "Fgf8", "Fgf10", "Uty")
vst.subset <- vst.c19[row.names(vst.c19) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D3 Target genes")


ddsHTSeq.mm.c19 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","N56.D3", "WT.D3")]
ddsHTSeq.mm.c19$treatment <- factor(as.character(ddsHTSeq.mm.c19$treatment))
ddsHTSeq.mm.c19$treatment <- relevel(ddsHTSeq.mm.c19$treatment, ref = "WT.D3")
dds.c19 <- DESeq(ddsHTSeq.mm.c19)

res.c19 <- results(dds.c19, filterFun=ihw, lfcThreshold = 1)

res.sig.c19 <- res.c19[!is.na(res.c19$padj),]
res.sig.c19 <- res.sig.c19[abs(res.sig.c19$log2FoldChange) > 1,]
res.sig.c19 <- res.sig.c19[res.sig.c19$padj < 0.1,]
res.sig.c19 <- res.sig.c19[order(res.sig.c19$pvalue),]
res.c19.df <- data.frame(res.c19)
res.sig.c19.df <- data.frame(res.sig.c19)
res.sig.c19.replicates <- c("WT_D3_1", "WT_D3_2", "WT_D3_3", "WT_D3_4", "N114_D3_1", "N114_D3_2", "N114_D3_3", "N114_D3_4", "N56_D3_1", "N56_D3_2", "N56_D3_3", "N56_D3_4")
WriteXLS(c("res.c19.df","res.sig.c19.df", "res.sig.c19.replicates"),"D3_replicates_Threshold_1_RNAseq_mESC_Epcam_N114.D3_N56.D3_vs_WT.D3.xlsx", SheetNames = c("all.results","stat.sig.results","res.sig.c19.replicates"), row.names = T)

heatmap_matrix <- as.matrix(hmd)
write.xlsx(heatmap_matrix, file = "heatmap_values.xlsx")

file_name <- "heatmap_data.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "heatmap_data")
writeData(wb, "heatmap_data", hmd)
saveWorkbook(wb, file_name)


matrix_data <- as.matrix(hmd)
summary(hmd)
summary_data <- summary(hmd)
wb <- createWorkbook()
addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_data)
saveWorkbook(wb, file = "summary_hmd.xlsx")
write.csv(summary_data, file = "summary_hmd.csv", row.names = T)
write.table(summary_data, file = "summary_hmd.txt", row.names = T, sep = "\t")

matrix_subset <- subset(hmd, rownames(hmd) %in% genes_of_interest)
summary(matrix_subset)
summary(data.frame(hmd[genes_of_interest,]))

for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the summary of the subsetted matrix
  print(paste(gene,":", sep = " "))
  print(summary(matrix_subset))
}

for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Print the gene name and mean values
  print(paste(gene,":", sep = " "))
  print(mean_values)
}

mean_values_list <- list()
for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the mean values of the subsetted matrix
  mean_values <- colMeans(matrix_subset)

  # Add the gene name and mean values to a list
  mean_values_list <- rbind(mean_values_list, c(gene,mean_values))
}
# Write the list to a csv file
write.csv(mean_values_list, file = "mean_values.csv", row.names = FALSE)


for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)
  # Get the summary of the subsetted matrix
  summary_data <- paste(gene,":", summary(matrix_subset), sep = " ")
  # Write the data to a tab-separated values file
  write.csv(summary_data, file = paste(gene,".csv", sep = ""))
}


# Create an empty data frame to store the summary of each gene
gene_summary <- data.frame()

# Loop through each gene
for (gene in genes_of_interest){
  # Subset the matrix to only include the gene of interest
  matrix_subset <- subset(hmd, rownames(hmd) == gene)

  # Get the summary of the subsetted matrix
  gene_summary <- rbind(gene_summary, summary(matrix_subset))
}

# Write the data frame to an excel file
write.table(gene_summary, file = "gene_summary.xls", sep = "\t", row.names = F, col.names = T)
write.xlsx(gene_summary, file = "gene_summary.xlsx")


genes_of_interest <- c("Gata6", "Cdkn2a", "T", "Dkk1", "Sall3", "Foxa2", "Wnt11", "Mesp1", "Eomes", "Mixl1", "Fgf8", "Fgf10", "Uty")
vst.subset <- vst.c19[row.names(vst.c19) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "D3 Target genes")

# Select the treatment groups of interest
treatment_of_interest <- c("N56.D3", "N114.D3", "WT.D3")
dds.c19_subset <- dds.c19[ ,colData(dds.c19)$treatment %in% treatment_of_interest]
genes_of_interest <- c("Gata6", "Cdkn2a", "T", "Dkk1", "Sall3", "Foxa2", "Wnt11", "Mesp1", "Eomes", "Mixl1", "Fgf8", "Fgf10", "Uty")
vst.subset <- vst(dds.c19_subset)[row.names(vst(dds.c19_subset)) %in% genes_of_interest,]
hmd <- assay(vst.subset)
colnames(hmd) <- dds.c19_subset$treatment
pheatmap(hmd, col = bluered(50), cluster_rows = T, scale = "row", main = "Selected treatment groups and genes of interest")

# calculate mean of replicates for each gene and treatment
mean.values <- rowMeans(hmd)

# subset data to include only treatment groups of interest and genes of interest
treatment.of.interest <- c("N56.D3", "N114.D3", "WT.D3")
hmd.subset <- hmd[,colnames(hmd) %in% treatment.of.interest]
genes.of.interest <- c("Gata6", "Cdkn2a", "T", "Dkk1", "Sall3", "Foxa2", "Wnt11", "Mesp1", "Eomes", "Mixl1", "Fgf8", "Fgf10", "Uty")
hmd.subset <- hmd.subset[row.names(hmd) %in% genes.of.interest,]

# plot pheatmap with the subsetted data
pheatmap(hmd.subset, col = bluered(50), cluster_rows = T, scale = "row", main = "D3 Target genes")
