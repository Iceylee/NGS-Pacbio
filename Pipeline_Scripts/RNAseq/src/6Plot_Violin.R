#!/usr/bin/env Rscript
args<-commandArgs(T)

# test if there are 4 arguments: if not, return an error
if (length(args)!=4) {
  print ("4 arguments must be supplied:")
  print ("FPKM_file coldata_file outPath whether_ref")
  print ("RSEM.FPKM.matrix ../R_input/colData.csv ./ FALSE")
  stop(call.=FALSE)
} 


FPKM_file = args[1] #"RSEM.FPKM.matrix"
coldata_file=args[2] #"colData.csv"
Path = args[3] #"./"
whether_ref = args[4] #FALSE

###library
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))


###Load data
if (whether_ref){
    df <- read_tsv(FPKM_file)
}else{
    df <- read_csv(FPKM_file)
    colnames(df)[1]="GeneID"
}

colData=read_csv(coldata_file)
colnames(colData) = c("sample","condition")

#按分组算每组的FPKM均值
groups = as.character(unique(colData$condition))
getGroupMean <- function(group){
    
    rows= dplyr::filter(colData, condition == group) %>% select(sample) %>% unlist
    df[[group]] <<- rowMeans(df[rows], na.rm = TRUE)
    return (1)

    }
sapply(groups, getGroupMean)


#melt转换
df <- select(df,groups,GeneID)
df <- melt(df, id = "GeneID")
colnames(df)  <- c("GeneID", "Group", "FPKM")
#log
df <- mutate(df, logFPKM=log10(FPKM+1))


###plotting
options(scipen=200)
p <- ggplot(df, aes(x = Group, y = logFPKM, fill = Group)) + 
geom_violin(alpha = 0.5) +  scale_y_log10() +
geom_boxplot(width=0.1, aes(fill=Group)) +
theme_classic() +
ylab("log10(FPKM+1)") +
xlab(NULL)
#由于scale_y_log10()，作图去掉了FPKM为0的值

###output
ggsave(file=paste(Path,"/violin_plot.pdf",sep=""),p, width=6, height=6, units="in")

ggsave(file=paste(Path,"/violin_plot.png",sep=""),p, width=6, height=6, units="in")