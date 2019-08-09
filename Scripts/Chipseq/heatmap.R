#!/usr/bin/env Rscript

args<-commandArgs(T)

group <- args[1]

df <- read.table(paste(group,"_Enrichment.txt",sep=""),header=T,row.names=1)
library("pheatmap")
library("genefilter")
topVarGenes <- head(order(rowVars(df), decreasing = TRUE), 1000)
df  <- df[ topVarGenes, ]
df <- log10(df + 1)
# anno = data.frame(Condition = colData$condition)
# anno$Condition = as.factor(anno$Condition)
# rownames(anno) = rownames(colData)
pdf(paste(group,"_heatmap.pdf",sep=""))
pheatmap(df,cluster_cols = FALSE,show_rownames = F)
dev.off()
