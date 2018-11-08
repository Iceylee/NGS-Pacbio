library(dplyr)
library(DESeq2)

########################################################
#####################4.deseq2-diff-gene#################
########################################################
load(paste(path1,"1Data_Input.RData",sep="/"))
groups = paste(exp_group,"vs",base_group,sep="")

#差异分析结果 multiple use contrast
res = results(dds,contrast=c("condition",exp_group,base_group))  
res_df = data.frame("id"=rownames(res),res)

regSig2 <- subset(res, abs(log2FoldChange) > 1)
sorted_regSig2 = regSig2[with(regSig2, order(-log2FoldChange)), ]
sig = data.frame("id"=rownames(regSig2),regSig2)
filename = paste(groups,"_sig_genes_exprData.txt",sep="")
write.table(sig, file= paste(path2,filename,sep="/"), sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

##scatter plot
res_df$CK=res_df$baseMean * 2 / (1 + 2^res_df$log2FoldChange)
res_df$Treat= res_df$CK * (2^res_df$log2FoldChange)


color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
res_df$color <- ifelse(abs(res_df$log2FoldChange)>=1,ifelse(res_df$log2FoldChange > 1 ,'red','blue'),'gray')

p <- ggplot(res_df, aes(x = Treat, y = CK,col = color,alpha=0.8)) +
  geom_point(size=1) +
  scale_color_manual(values = color) +
  labs(x="Treat",y="CK")+
  theme(legend.position = "none",
        panel.grid=element_blank()) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14)) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10') #对数坐标轴




p

filename=paste("lncRNA2_scatter_plot.pdf",sep="")
ggsave(filename,p, width=8, height=6, units="in")
