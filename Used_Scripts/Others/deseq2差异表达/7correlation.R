setwd("~/work/11月/1212GO&KEGG/")
df <- read.table(file = "norm-matrix-deseq2.txt",sep = '\t',header = T,stringsAsFactors=F,row.names = 1)
cor.df <- cor(df[,unlist(lapply(df, is.numeric))])

library(pheatmap)
pdf(file="heatmap_cor.pdf")
pheatmap(cor.df,legend = TRUE,display_numbers = TRUE, number_format = "%.2f", number_color="black")
dev.off()

#####RPKM
library(reshape2)
rpkm_df <- read.table(file = "AllSamplesRPKMValue.txt",sep = '\t',header =F,stringsAsFactors=F)
dff <- rpkm_df[,1:3]
rpkm_tab <- dcast(dff, V2 ~ V1)
row.names(rpkm_tab) <- rpkm_tab$V2
df2 <- rpkm_tab[,-1]

cor.df2 <- cor(df2[,unlist(lapply(df2, is.numeric))])
library(pheatmap)
pdf(file="heatmap_cor.pdf")
pheatmap(cor.df2,legend = TRUE,display_numbers = TRUE, number_format = "%.2f", number_color="black")
dev.off()
