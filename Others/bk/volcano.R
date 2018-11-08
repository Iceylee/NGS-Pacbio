log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("dplyr")

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