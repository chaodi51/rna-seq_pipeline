## Setup contrast as Coculture vs Cdiff or Efaecalis, oupout tables with the standard DESeq2 result format including " 
## baseMean log2FoldChange lfcSE stat pvalue padj" plus the “Log fold change shrinked” normalized readcounts 

library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["count_table"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@params[["sample_contrast"]], header=TRUE, check.names=FALSE) 
# cts <- read.table("../results/tables/all_readCount.tsv", header=TRUE,row.names="gene", check.names=FALSE)
# coldata <- read.table("sample_contrast.tsv",header=TRUE, check.names=FALSE)
rownames(coldata)=coldata$sample

dds <- DESeqDataSetFromMatrix(countData=cts, colData=coldata, design = ~ condition)
dds$condition <- relevel(dds$condition, "WT") # use "WT" as the reference

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 10, ]
# normalization and pre-processing
dds <- DESeq(dds, parallel=parallel)

# raw count normalization
norm_counts <- counts(dds, normalized=TRUE) 
# count transformation, log2 scale, either rlog or vst
vsd <- vst(dds, blind=FALSE)

# get the current contrast/coculture_vs_mono from snakemake output, e.g., "Coculture_vs_Cdiff"
output_file <- snakemake@output[["table"]]
# output_file <- "LPS_vs_WT.diffexp.tsv"
comp = gsub(".diffexp.tsv", "", tail(unlist(strsplit(output_file, "/")),1))

res <- results(dds, contrast = c("condition", unlist(strsplit(comp, "_vs_"))[1], unlist(strsplit(comp, "_vs_"))[2]), parallel = parallel)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast = c("condition", unlist(strsplit(comp, "_vs_"))[1], unlist(strsplit(comp, "_vs_"))[2]), res=res, type="ashr")
# report # up/down genes
paste0("padj<=0.01 & log2FoldChange>=1: # Up = ", length(res[which(res$padj<=1e-2 & res$log2FoldChange>1),]$padj),
       "  # Down = ", length(res[which(res$padj<=1e-2 & res$log2FoldChange < -1),]$padj))

# extract the current cell_type samples
df_vsd = as.data.frame(assay(vsd))
# merge with normalized count data and output the table
resdata <- merge(df_vsd, as.data.frame(res), by="row.names",sort=FALSE)
names(resdata)[1] <- "Gene"

write.table(resdata, file=snakemake@output[["table"]], sep="\t", quote=FALSE, row.names=FALSE)
# write.table(resdata, file="../results/diffexp/LPS_vs_WT.diffexp.tsv", sep="\t", quote=FALSE, row.names=FALSE)
## basic plots for Data quality assessment
# M-A plot, points are red with padj < 0.1, points fall out of the window are open triangles 
pdf(snakemake@output[["ma_plot"]])
# pdf("ma_plot.pdf")
plotMA(res, main=comp, colLine="red")
dev.off()

# volcano plot using EnhancedValcano
pdf(snakemake@output[["volcano_plot"]])
# pdf("volcano-plot.pdf",10,10)
EnhancedVolcano(resdata,
                lab = resdata$Gene,
                x = 'log2FoldChange',
                y = 'padj',
                title="LPS vs WT",
                subtitle = "",
                pCutoff = 1e-2,
                FCcutoff = 1,
                pointSize = 2,
                labSize = 4) + 
                theme(plot.title = element_text(hjust = 0.5)) 

dev.off()
