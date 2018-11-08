library(dplyr)
library(pheatmap)
setwd("~/work/11月/1121cuffdiff/")

#1. deseq2 raw count 
###dds 未标准化
load(file = "1Data_Input.RData")

rld <- rlog(dds, blind = FALSE)
#######前20样本差异较大的，做聚类热图。不同样本中的表达情况
library("genefilter") 
library("pheatmap")

topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30) 
mat  <- assay(rld)[ topVarGenes, ] 
mat  <- mat - rowMeans(mat) 

#从coldata导入condition信息
# coldata_file = "colData.csv"
# colData = read.csv(coldata_file, header=T, sep="\t", row.names=1 )
# anno = data.frame(Condition = colData$condition)

#手动写condition信息
anno = data.frame(Condition = c("BZ", "BZ", "V199", "V199"))

anno$Condition = as.factor(anno$Condition)

rownames(anno) = rownames(colData)

pdf(file="heatmap_top30.pdf")
#cluster_cols=F
pheatmap(mat,annotation_col = anno)
dev.off()



###2.fpkm做热图

fpkm_data <- read.table(file = "Cuffdiff/genes.fpkm_tracking", header=TRUE, sep="\t", stringsAsFactors = F)

fpkm_data2 <- select(fpkm_data,Pal_12dpi_FPKM,Pal.103_12dpi_FPKM,gene_short_name)
colnames(fpkm_data2)[3] <- "gene" 


diff_data <- read.table(file = "gene_exp.diff", header=TRUE, sep="\t", stringsAsFactors = F)
sig_diff <- filter(diff_data,abs(log2.fold_change.) > 1,q_value < 0.05)
write.table(sig_diff, file="sig-genes-cuffdiff.txt", sep="\t",  row.name=F, col.names=TRUE,quote=FALSE)

result <- dplyr::left_join(sig_diff,fpkm_data2,by = "gene")

#get q value top 50
result_q <- arrange(result,q_value) 
result_q_top <- result_q[1:30,15:16] 
rownames(result_q_top) <- result_q[1:30,3]

#log变换
result_q_top_log <- log10(result_q_top + 1)


#添加了legend的标题 log10(FPKM+1)
pdf(file="heatmap_top30.pdf",onefile=FALSE))
pheatmap(result_q_top_log,scale="column", legend_breaks = c(0, 1, 2,max(result_q_top_log)),legend_labels = c("0", "1", "2","log10(FPKM+1)\n"),
         legend = TRUE)
dev.off()


##
# test <- matrix(rexp(200, rate=.1), ncol=20)
# colnames(test) = paste("Room", 1:20, sep = "")
# rownames(test) = paste("Building", 1:10, sep = "")


# pheatmap(test, legend_breaks = c(10, 20, 30, 40, max(test)), 
#          main = "", legend_labels = c("10", "20", "30", "40", "title\n"),
#          legend = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)
