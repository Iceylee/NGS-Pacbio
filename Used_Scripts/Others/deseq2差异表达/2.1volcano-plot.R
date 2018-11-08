########volcano2
library(ggplot2)
library(dplyr)

setwd("~/work/11月/1121cuffdiff/")

#1.cuffdiff数据输入

diff_data <- read.table(file = "gene_exp.diff", header=TRUE, sep="\t", stringsAsFactors = F)

#这里是cuffdiff的输出文件，包含所有的gene，第一步根据teststat筛选出不是na的基因
sig_df <- filter(diff_data,!is.na(test_stat))
colnames(sig_df)[10]<- "log2FoldChange"
colnames(sig_df)[13] <- "padj"

##2.deseq数据输入
load(file = "1Data_Input.RData")
#res = results(dds,contrast=c("condition","MS","V199"))  #多组的话 两两比对
df = data.frame("id"=rownames(res),res) #使用所有的基因，而不是筛选过的显著差异基因
sig_df <- filter(df,!is.na(padj)) #去除qadj为NA的数据

##3.作图

#adj<0.05 AND log2foldchange>1
#scale_color_manual 
color<- c(red = "red", gray = "gray", blue ="blue")
#add color column ,and condition
sig_df$color <- ifelse(sig_df$padj < 0.05 & abs(sig_df$log2FoldChange) >=1,ifelse(sig_df$log2FoldChange > 1 ,'red','blue'),'gray')

library(ggplot2)
p <- ggplot(sig_df, aes(x = log2FoldChange, y = -log10(padj),col = color)) +
  geom_point() +
  scale_color_manual(values = color) +
  labs(x="log2 (fold change)",y="-log10 (padj)")+
  geom_hline(yintercept = -log10(0.05), lty=4,col="grey",lwd=0.6) +
  geom_vline(xintercept = c(-1, 1), lty=4,col="grey",lwd=0.6) +
  theme(legend.position = "none",
        panel.grid=element_blank())
  #       +theme（axis.title = element_text(size = 16),
  #       axis.text = element_text(size = 14))+
  # scale_x_continuous(limits = c(-10, 10)) +
  # scale_y_continuous(limits = c(0,3)))

pdf(file="volcano.pdf")
p
dev.off()