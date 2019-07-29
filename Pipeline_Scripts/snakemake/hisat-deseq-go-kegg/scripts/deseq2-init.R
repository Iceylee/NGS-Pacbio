log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# Read in table
countData = read.csv(snakemake@input[["countTab"]], header=T, row.names=1)
colData = read.csv(snakemake@input[["colData"]], header=T,row.names = 1)
colData$condition = as.factor(colData$condition)
if (all(rownames(colData) %in% colnames(countData)) == F) {
  print ("rownames(colData) != colnames(countData)")
  print ("Please check.")
  stop(call.=FALSE)
}

# Create DESEq2 dataset.
dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition) #~在R里面用于构建公式对象，~左边为因变量，右边为自变量。

#pre filter
dds = dds[ rowSums(counts(dds)) > 1 ,] 
#Set the reference to be compared
#dds$condition <- relevel(dds$condition,snakemake@params[["base_group"]])

# Run deseq
dds <- DESeq(dds) #标准化
# Get normalized counts and write this to a file
nc = counts(dds,normalized=TRUE)
# Turn it into a dataframe to have proper column names.
dt = data.frame("id"=rownames(nc),nc)
colnames(dt) <- c("id",colnames(nc))
# Save the normalize data matrix.
write.table(dt, file=snakemake@output[["countMatrix"]], sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)

saveRDS(dds, file=snakemake@output[["rds"]])
