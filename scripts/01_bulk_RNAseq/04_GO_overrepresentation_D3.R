## =============================================================================
## Related : Gong N, Gouda M, et al. Cell Death & Disease 17, 389 (2026),
##           https://doi.org/10.1038/s41419-026-08734-w (see README: Disclaimer)
## Module  : 01 | Bulk 3'-RNA-seq of WT and Epcam-/- mESC (clones #56, #114), D0-D10 (GEO: GSE293121)
## Script  : 04_GO_overrepresentation_D3.R
## Purpose : GO biological-process over-representation (clusterProfiler::enrichGO) of genes
##           up-regulated in pooled Epcam-/- (#56/#114) vs WT EBs at D3 (DESeq2 + IHW,
##           lfcThreshold = 1; padj < 0.05, baseMean > 50, log2FC > 0.5).
## Input   : data/bulk_RNAseq/featurecounts_mESC_EpCAM_bulkRNAseq.Rdata
## Output  : results/01_bulk_RNAseq/GO_BP_D3_up_KO_vs_WT_dotplot.png
## =============================================================================

## ---- Packages -----------------------------------------------------------------

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(Rsubread)
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
library(gdata)
library(knitr)
library(gtools)

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


## ---- Pooled KO (#56 + #114) vs WT at D3 (c19) ------------------------------
ddsHTSeq.mm.c19 <- ddsHTSeq.mm[,ddsHTSeq.mm$treatment%in%c("N114.D3","N56.D3", "WT.D3")]
ddsHTSeq.mm.c19$treatment <- factor(as.character(ddsHTSeq.mm.c19$treatment))
ddsHTSeq.mm.c19$treatment <- relevel(ddsHTSeq.mm.c19$treatment, ref = "WT.D3")
dds.c19 <- DESeq(ddsHTSeq.mm.c19)

res.c19 <- results(dds.c19, filterFun=ihw, lfcThreshold = 1)

summary(res.c19)

res.c19.plot <-  res.c19[!is.na(res.c19$padj),]

v.c19 <- data.frame(res.c19.plot[,c(2,6)])
colnames(v.c19)[2] <- "FDR"
ev.c19 <- EnhancedVolcano(v.c19,
                          lab = rownames(v.c19),
                          x = 'log2FoldChange',
                          y = 'FDR',
                          title = 'N114 D3 N56 D3 versus WT D3',
                          pCutoff = 0.25,
                          FCcutoff = 0,
                          pointSize = 3.0,
                          labSize = 3.0,
                          ylab=bquote(~-Log[10] ~ italic(FDR)))
ev.c19
res.sig.c19 <- res.c19[!is.na(res.c19$padj),]
res.sig.c19 <- res.sig.c19[abs(res.sig.c19$log2FoldChange) > 1,]
res.sig.c19 <- res.sig.c19[res.sig.c19$padj < 0.1,]
res.sig.c19 <- res.sig.c19[order(res.sig.c19$pvalue),]
res.c19.df <- data.frame(res.c19)
res.sig.c19.df <- data.frame(res.sig.c19)

res.sig.c19


## ---- GO:BP over-representation of up-regulated D3 genes ------------------
sigs <- na.omit(res.sig.c19)
sigs <- sigs [sigs$padj < 0.05 & sigs$baseMean > 50,]
sigs
genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5,])
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP")
as.data.frame(GO_results)

fit <- plot(dotplot(GO_results, showCategory = 50))
fit

fit <- plot(barplot(GO_results, showCategory = 50))
fit

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
as.data.frame(GO_results)
fit <- plot(dotplot(GO_results, showCategory = 20))
dotplot(GO_results)

png("GO_BP_D3_up_KO_vs_WT_dotplot.png", res = 250, width = 1200, height = 1000)
print(fit)
dev.off()
