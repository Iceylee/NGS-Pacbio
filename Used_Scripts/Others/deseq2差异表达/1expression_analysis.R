setwd("~/work/11月/1212GO&KEGG/")

# load count table and coldata
count_table = "CountMatrix4DESeq.csv"
coldata_file = "colData.csv"


# Load the library.
library(DESeq2)
library(dplyr)


#sep NEED CHANGE
###1.准备
countData = read.table(count_table, header=TRUE, sep=",", row.names=1)
View(countData)

colData = read.csv(coldata_file, header=T, sep="\t", row.names=1 )
View(colData)

colData$condition = as.factor(colData$condition)
colnames(countData) <- rownames(colData)
#check
all(rownames(colData) %in% colnames(countData))


###2.Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
#pre filter
dds = dds[ rowSums(counts(dds)) > 1 ,] 
#Set the reference to be compared
####NEED CHANGE#######
dds$condition = relevel(dds$condition,"MS")


###3.Run deseq
dds = DESeq(dds) #标准化
#提取counts数据
nc = counts(dds,normalized=TRUE)
#转df，加行名
dt = data.frame("id"=rownames(nc),nc)
#输出标准化count矩阵
write.table(dt, file="norm-matrix-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

###4.差异分析结果 
##4.1 多个实验组 定义参数contrast
#res = results(dds,contrast=c("condition","MS","BZ")) 
##4.2 仅一个实验组和一个对照组 
res = results(dds）
# Sort the results data frame by the padj and foldChange columns.
sorted = res[with(res, order(padj, -log2FoldChange)), ]
#转df，添加行名
sorted.df = data.frame("id"=rownames(sorted),sorted)
#输出全部基因差异分析结果
write.table(sorted.df, file="deseq2.txt", row.names = FALSE,sep="\t", quote=FALSE)

###5.筛选显著差异基因
#padj（<0.05)和log2 fold (>1)
regSig <- subset(res, padj < 0.05)
regSig2 <- subset(regSig, abs(log2FoldChange) > 1)
sorted_regSig2 = regSig2[with(regSig2, order(-log2FoldChange)), ]
sig = data.frame("id"=rownames(regSig2),regSig2)
#输出显著差异基因分析结果
write.table(sig, file="sig-genes-deseq2.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)


save(dds, res, sorted.df, dt, sig, file="1Data_Input.RData")


