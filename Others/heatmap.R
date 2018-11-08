log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

dds <- readRDS(snakemake@input[["rds"]])
########################################################
#heatmap
########################################################
rld <- rlog(dds, blind = FALSE)
#######前20样本差异较大的，做聚类热图。不同样本中的表达情况
library("genefilter") 
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30)
mat  <- assay(rld)[ topVarGenes, ] 
mat  <- mat - rowMeans(mat) 
anno = data.frame(Condition = colData$condition)
anno$Condition = as.factor(anno$Condition)
rownames(anno) = rownames(colData)
pdf(snakemake@output[[1]],onefile=FALSE)
#cluster_cols = F
pheatmap(mat,annotation_col = anno,cluster_cols = FALSE)
dev.off()