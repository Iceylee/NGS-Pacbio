load(file = "1Data_Input.RData")

rld <- rlog(dds, blind = FALSE)
#######前20样本差异较大的，做聚类热图。不同样本中的表达情况
library("genefilter") 
library("pheatmap")

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 20) 
mat  <- assay(rld)[ topVarGenes, ] 
mat  <- mat - rowMeans(mat) 

anno = data.frame(Condition = c("BZ", "BZ", "V199", "V199"))
anno$Condition = as.factor(anno$Condition)

rownames(anno) = rownames(colData)

pdf(file="heatmap_top20.pdf")
pheatmap(mat,annotation_col = anno)
dev.off()

#######火山图
# sig <- res[which(res$pvalue < 0.05),] 
# cols <- rep("#000000", nrow(res)) 
# cols[res$log2FoldChange >= 2] <- "#CC0000" 
# cols[res$log2FoldChange <= -2] <- "#0000CC" 
# 
# pdf(file="volcano.pdf")
# plot(sig$log2FoldChange, -log10(sig$pvalue), pch=16, cex=0.75,col=cols,las=1 )
# dev.off()

########volcano2
library(ggplot2)
sig <- res[which(res$padj < 0.05),] 
sig_df <- as.data.frame(sig)

#scale_color_manual 
color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
sig_df$color <- ifelse(sig_df$padj < 0.05 & abs(sig_df$log2FoldChange) >= 1,ifelse(sig_df$log2FoldChange > 1 ,'red','blue'),'gray')

p <- ggplot(sig_df, aes(x = log2FoldChange, y = -log10(padj),col = color)) +
  geom_point() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (padj)")+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14))
pdf(file="volcano.pdf")
p
dev.off()
