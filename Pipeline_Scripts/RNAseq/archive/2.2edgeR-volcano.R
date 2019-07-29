###下载安装edgeR包
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
library("edgeR")
library('ggplot2')

########################################################
#####################1.edgeR-diff-gene#################
########################################################

###读取数据
countData = read.table(count_table, header=TRUE, sep=",", row.names=1)
colData = read.csv(coldata_file, header=T,row.names = 1)
groups = paste(exp_group,"vs",base_group,sep="")

#提取
colData$smp = rownames(colData)
base_smp = colData[colData$condition==base_group,]$smp
exp_smp = colData[colData$condition==exp_group,]$smp

#进行分组
rawdata <- countData[,c(base_smp,exp_smp)] #base should be in first column
group <- factor(c(base_smp,exp_smp))

###过滤与标准化
y <- DGEList(counts=rawdata,genes=rownames(rawdata),group = group)
###TMM标准化
y<-calcNormFactors(y)
y$samples
###推测离散度,根据经验设置，若样本是人，设置bcv = 0.4，模式生物设置0.1.
#bcv <- 0.1
bcv <- 0.2
#bcv <- 0.4
et <- exactTest(y, dispersion=bcv^2)
topTags(et)
summary(de <- decideTestsDGE(et))

###导出数据
DE <- et$table
DE$significant <- as.factor(DE$PValue<0.05 & abs(DE$logFC) >1)
#write.table(DE,file="edgeR_all2",sep="\t",na="NA",quote=FALSE)

filename = paste(groups,"_all_genes_exprData.txt",sep="")
write.table(DE, file= paste(path2,filename,sep="/"), sep="\t",  row.name=TRUE, col.names=TRUE,quote=FALSE)

filename = paste(groups,"_sig_genes_exprData.txt",sep="")
DE_sig <- DE[DE$significant=="TRUE",]
write.table(DE_sig, file= paste(path2,filename,sep="/"), sep="\t",  row.name=TRUE, col.names=TRUE,quote=FALSE)

########################################################
#####################2.MA plot##########################
########################################################
filename = paste(groups,"_MA_plot.pdf",sep="")
pdf(file = paste(path2,filename,sep="/"))
detags <- rownames(y)[as.logical(DE$significant)]
#detags <- rownames(y)[as.logical(de)];
plotSmear(et, de.tags=detags)
abline(h=c(-1, 1), col="blue");
dev.off()

########################################################
#####################3.volcano plot#####################
########################################################

df <- DE
sig_df <- df
#df = data.frame("id"=rownames(res),res) #使用所有的基因，而不是筛选过的显著差异基因
#sig_df <- filter(df,!is.na(padj)) #去除qadj为NA的数据
names(sig_df) <- c("log2FoldChange","logCPM","padj","significant")

#adj<0.05 AND log2foldchange>1
#scale_color_manual 
color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
sig_df$color <- ifelse(sig_df$padj < 0.05 & abs(sig_df$log2FoldChange) >=1,ifelse(sig_df$log2FoldChange > 1 ,'red','blue'),'gray')

library(ggplot2)
p2 <- ggplot(sig_df, aes(x = log2FoldChange, y = -log10(padj),col = color)) +
  geom_point() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (padj)")+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank())

filename = paste(groups,"_volcano_plot.pdf",sep="")
ggsave(file=paste(path2,filename,sep="/"),p2, width=6, height=6, units="in")


