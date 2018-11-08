library(dplyr)
library(DESeq2)

########################################################
#####################4.deseq2-diff-gene#################
########################################################
load(paste(path1,"1Data_Input.RData",sep="/"))
groups = paste(exp_group,"vs",base_group,sep="")

#差异分析结果 multiple use contrast
res = results(dds,contrast=c("condition",exp_group,base_group))  
# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]
# Turn it into a dataframe to have proper column names.
sorted.df = data.frame("id"=rownames(sorted),sorted)
# Write the table out.
filename=paste(groups,"_all_genes_exprData.txt",sep="")
write.table(sorted.df, file=paste(path2,filename,sep="/"), row.names = FALSE,sep="\t", quote=FALSE)

########significantly different genes
#padj（<0.05)和log2 fold (>1)
regSig <- subset(res, padj < 0.05)
regSig2 <- subset(regSig, abs(log2FoldChange) > 1)
sorted_regSig2 = regSig2[with(regSig2, order(-log2FoldChange)), ]
sig = data.frame("id"=rownames(regSig2),regSig2)
filename = paste(groups,"_sig_genes_exprData.txt",sep="")
write.table(sig, file= paste(path2,filename,sep="/"), sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

########################################################
#####################5.volcano##########################
########################################################
df = data.frame("id"=rownames(res),res)
sig_df <- filter(df,!is.na(padj))
#adj<0.05 AND log2foldchange>1
#scale_color_manual 
color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
sig_df$color <- ifelse(sig_df$padj < 0.05 & abs(sig_df$log2FoldChange) >=1,ifelse(sig_df$log2FoldChange >= 1 ,'red','blue'),'gray')
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

filename=paste(groups,"_volcano_plot.pdf",sep="")
ggsave(file=paste(path2,filename,sep="/"),p, width=6, height=6, units="in")



