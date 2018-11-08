log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


dds <- readRDS(snakemake@input[["rds"]])
########################################################
#deseq2-diff-gene
########################################################

#差异分析结果 multiple use contrast
contrast <- c("condition", snakemake@params[["contrast"]])
res = results(dds,contrast=contrast)  
# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]
# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)
# Write the table out.
write.table(sorted.df, file=snakemake@output[["all_tab"]], row.names = FALSE,sep="\t", quote=FALSE)

########significantly different genes
#padj（<0.05)和log2 fold (>1)
regSig <- subset(res, padj < 0.05)
regSig2 <- subset(regSig, abs(log2FoldChange) > 1)
sorted_regSig2 = regSig2[with(regSig2, order(-log2FoldChange)), ]
sig = data.frame("id"=rownames(regSig2),regSig2)
write.table(sig, file= snakemake@output[["sig_tab"]], sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)