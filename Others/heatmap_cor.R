log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

########################################################
#correlation
########################################################
df <- read.table(file = snakemake@input[[1]],sep = '\t',header = T,stringsAsFactors=F,row.names = 1)
colnames(df) <- rownames(colData)
cor.df <- cor(df[,unlist(lapply(df, is.numeric))])

library(pheatmap)
pdf(snakemake@output[[1]],onefile=FALSE)
pheatmap(cor.df,legend = TRUE,display_numbers = TRUE, number_format = "%.2f", number_color="black",cluster_cols = FALSE)
dev.off()

