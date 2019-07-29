log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("dplyr")

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

########################################################
#####################5.volcano##########################
########################################################

df = data.frame("id"=rownames(res),res)
sig_df <- filter(df,!is.na(padj))
#adj<0.05 AND log2foldchange>1
#scale_color_manual 
color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
sig_df$color <- ifelse(sig_df$padj < 0.05 & abs(sig_df$log2FoldChange) >=1,ifelse(sig_df$log2FoldChange > 1 ,'red','blue'),'gray')
##
library(ggplot2)
rm(p)
p <- ggplot(sig_df, aes(x = log2FoldChange, y = -log10(padj),col = color)) +
  geom_point() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (padj)")+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank()) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14))

ggsave(file=snakemake@output[["plot"]],p, width=6, height=6, units="in")
