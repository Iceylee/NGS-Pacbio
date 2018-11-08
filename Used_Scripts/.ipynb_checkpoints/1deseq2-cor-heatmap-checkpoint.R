

########################################################
#####################1.deseq2-normalize#################
########################################################
# Load the library.
library(DESeq2)
library(dplyr)

countData = read.table(count_table, header=TRUE, sep=",", row.names=2)
countData = countData[,-1] #delete first column(entrezID)
colData = read.csv(coldata_file, header=T)
colData$condition = as.factor(colData$condition)
colnames(countData) <- colData[,1]
#all(rownames(colData) %in% colnames(countData))

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。

#pre filter
dds = dds[ rowSums(counts(dds)) > 1 ,] 
#Set the reference to be compared
dds$condition <- relevel(dds$condition,base_group)

# Run deseq
dds <- DESeq(dds) #标准化
# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)
# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)
colnames(dt) <- c("id",colnames(nc))
# Save the normalize data matrix.
write.table(dt, file=paste(path1,"norm-count-matrix.txt",sep="/"), sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

save(dds, file=paste(path1,"1Data_Input.RData",sep="/"))

########################################################
#####################2.correlation######################
########################################################
df <- read.table(file = paste(path1,"norm-count-matrix.txt",sep="/"),sep = '\t',header = T,stringsAsFactors=F,row.names = 1)
colnames(df) <- colData[,1]
cor.df <- cor(df[,unlist(lapply(df, is.numeric))])

library(pheatmap)
pdf(file=paste(path1,"heatmap_cor.pdf",sep="/"),onefile=FALSE)
pheatmap(cor.df,legend = TRUE,display_numbers = TRUE, number_format = "%.2f", number_color="black")
dev.off()

########################################################
#####################3.heatmap##########################
########################################################
rld <- rlog(dds, blind = FALSE)
#######前20样本差异较大的，做聚类热图。不同样本中的表达情况
library("genefilter") 
library("pheatmap")
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 30)
mat  <- assay(rld)[ topVarGenes, ] 
mat  <- mat - rowMeans(mat) 
anno = data.frame(Condition = colData$condition)
anno$Condition = as.factor(anno$Condition)
rownames(anno) = colData$Sample
pdf(file=paste(path1,"heatmap_top30.pdf",sep="/"),onefile=FALSE)
#cluster_cols = F
pheatmap(mat,annotation_col = anno)
dev.off()



