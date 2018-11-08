
suppressMessages(library("dplyr"))
suppressMessages(library('ggplot2'))

########################################################
#####################1.Gfold-diff-gene#################
########################################################

###读取数据
countData = read.table(count_table, header=TRUE, sep=",", row.names=1)
colData = read.csv(coldata_file, header=T,row.names = 1)
countData <- countData[,rownames(colData)]
groups = paste(exp_group,"vs",base_group,sep="")

#将列名：编号转成样本名称
colnames(countData) = colData$condition


#将第二，四，五列补齐为1.
geneSum = dim(countData)[1]
placeHolder = rep(1,geneSum)

#创建实验组和对照组的cnt文件。
cntData_exp = data.frame(rownames(countData),placeHolder, select(countData,!!sym(exp_group)), placeHolder, placeHolder)
write.table(cntData_exp, file=paste(path2,"/",exp_group,".read_cnt",sep=""), col.names = FALSE,row.names = FALSE, sep="\t", quote=FALSE)

cntData_base = data.frame(rownames(countData),placeHolder, select(countData,!!sym(base_group)), placeHolder, placeHolder)
write.table(cntData_base, file=paste(path2,"/",base_group,".read_cnt",sep=""), col.names = FALSE,row.names = FALSE, sep="\t", quote=FALSE)

#gfold
diffFileName = paste(path2,"/",exp_group,"vs",base_group,".diff",sep="")
system(sprintf("gfold diff -s1 %s -s2 %s -suf .read_cnt -sc 0.05 -o %s",paste(path2,"/",base_group,sep=""),paste(path2,"/",exp_group,sep=""),diffFileName),inter = TRUE)

diffData = read.table(diffFileName, header=FALSE, sep="\t", row.names=1,comment.char = "#")

#calculate baseMean
countMean = rowMeans(countData[c(exp_group, base_group)])
#exact gfold value
allGeneData = data.frame(Gene_ID = rownames(diffData), baseMean = countMean, GFOLD = diffData$V3, log2FoldChange = diffData$V3)
sigGeneData = filter(allGeneData,abs(GFOLD)>=1)

#output
write.table(allGeneData, file=paste(path2,"/",groups,"_all_genes_exprData.txt",sep=""), row.names = FALSE, sep="\t", quote=FALSE)
write.table(sigGeneData, file=paste(path2,"/",groups,"_sig_genes_exprData.txt",sep=""), row.names = FALSE, sep="\t", quote=FALSE)

########################################################
#####################2.MA-plot##########################
########################################################

allGeneData$logbaseMean = log2(allGeneData$baseMean+1)

#add color column ,and condition
color<- c(red = "red", gray = "gray")
allGeneData$color <- ifelse(abs(allGeneData$GFOLD) >=1,"red","gray")

p <- ggplot(allGeneData, aes(y = GFOLD, x = logbaseMean, col = color)) +
    geom_point() +scale_color_manual(values = color) +theme_light()+
    labs(y = "log2 Fold Change", x = "Mean Expression") +
    theme(legend.position = "none",panel.grid=element_blank()) +
    geom_hline(yintercept = c(-1, 1), lty=2,col="orange",lwd=0.6)+
    theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14)) +
    labs(title = paste(groups," MA plot",sep=""))+
  	theme(plot.title = element_text(hjust = 0.5))

ggsave(file=paste(path2,"/",groups,"_MA.pdf",sep=""),p, width=6, height=6, units="in")

ggsave(file=paste(path2,"/",groups,"_MA.png",sep=""),p, width=6, height=6, units="in")
