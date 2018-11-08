library(dplyr)
library(pheatmap)
setwd("~/work/11月/1121cuffdiff/")

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
result_q_top_log <- log10(result_q_top + 1)

pdf(file="heatmap_top30.pdf")
pheatmap(result_q_top_log, legend_breaks = c(0, 1, 2,max(result_q_top_log)),legend_labels = c("0", "1", "2","log10(FPKM+1)\n"),
         legend = TRUE)
dev.off()


##调参数用 legend标题可修改
test <- matrix(rexp(200, rate=.1), ncol=20)
colnames(test) = paste("Room", 1:20, sep = "")
rownames(test) = paste("Building", 1:10, sep = "")

pheatmap(test, legend_breaks = c(10, 20, 30, 40, max(test)), scale="column"，
         main = "", legend_labels = c("10", "20", "30", "40", "title\n"),
         legend = TRUE, cluster_rows = FALSE, cluster_cols = FALSE)

##
data<- read.table("221.txt",head = T,row.names=1,sep="\t")
pheatmap(data,
         scale="column",#按行均一化，"row","column" or "none"默认是"none"
         #clustering_distance_rows = "correlation",#聚类线长度优化
#        treeheight_row=40,#按行聚类树高
#        treeheight_col=40,#按列聚类树高
         cluster_cols=FALSE,#是否按列聚类
         cluster_rows=F,#是否按行聚类
         display_numbers=F,#是否在每一格上显示数据
         number_format="%.2f",#显示数据的格式，几位小数，或"\%.1e"，颜色number_color,大小fontsize_number
         fontsize_row=10,#行名称字体大小
         fontsize_col=15,#列名称字体大小
#格子大小
         cellwidth = 50,
         cellheight= 14,
         #main="ABC",#标题名称
         #gaps_row = c(10, 15),#插入缝隙，不能聚类！
         #cutree_row = 7,#按聚类分割
         show_colnames=TRUE,#是否显示列名，同理show_rownames
#定义颜色"navy", "white", "firebrick3"
         #color = colorRampPalette(c("green","white","red"),bias=2.5)(256),         
         color = colorRampPalette(c("MediumBlue","white","red"))(256),
         #border_color = "black",
#格子框颜色
         #legend = FALSE,#是否显示图例
         #legend_breaks = -1:4,#图例范围
         #filename = "test.pdf",#保存文件命名
)
