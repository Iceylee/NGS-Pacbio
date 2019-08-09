library("dplyr")

setwd("~/work/2018/0508差异基因lncRNA-2/")
countMat = read.table(file="norm-count-matrix.txt", header=T,stringsAsFactors = F)
countMat[,c(2,3,4,5,6,7)] = sapply(countMat[,c(2,3,4,5,6,7)],as.numeric)

countMat[countMat==0] = 1



countMat$Treat <- rowMeans(subset(countMat, select = c(2,3,4)), na.rm = TRUE)
countMat$CK = rowMeans(subset(countMat, select = c(5,6,7)), na.rm = TRUE)
countMat$log10Treat = log10(countMat$Treat)
countMat$log10CK = log10(countMat$CK)

countMat$log2FoldChange=log2(countMat$Treat/countMat$CK)

color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
countMat$color <- ifelse(abs(countMat$log2FoldChange)>=1,ifelse(countMat$log2FoldChange > 1 ,'red','blue'),'gray')

p <- ggplot(countMat, aes(x = Treat, y = CK,col = color,alpha=0.8)) +
  geom_point(size=1) +
  scale_color_manual(values = color) +
  labs(x="Treat",y="CK")+
  theme(legend.position = "none",
        panel.grid=element_blank()) +
  theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14)) +
  scale_y_continuous(trans = 'log10') +
  scale_x_continuous(trans = 'log10')

filename=paste("lncRNA2_scatter_plot.pdf",sep="")
ggsave(filename,p, width=8, height=6, units="in")

#sig gene list
sig_list <- countMat[abs(countMat$log2FoldChange)>=1,c(1,2)]
all_gene <- read.table(file="TreatvsCK_all_genes_exprData.txt", header=T,stringsAsFactors = F)
gene_list = left_join(sig_list,all_gene,by="id")
gene_list$F.1.1_Count=NULL #删除某列

write.table(gene_list,file = "lncRNA2_TreatvsCK_sig_genes_exprData.txt",sep = "\t",row.names = F,quote = F)