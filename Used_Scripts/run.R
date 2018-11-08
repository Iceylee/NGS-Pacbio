#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=7) {
	print ("7 arguments must be supplied:")
	print ("matrixFile groupFile baseGroup dbName keggName keyType keggType")
	print ("CountMatrix4DESeq.csv colData.csv NC AH134 hsa ENSEMBL ENTREZID")
  stop(call.=FALSE)
} 


count_table <- args[1] #"CountMatrix4DESeq.csv"
coldata_file <- args[2] #"colData.csv"
base_group <- args[3] #"NC"
dbname <- args[4] #'AH134'
kegg_org <- args[5]  #"hsa"
GO_KEY <- args[6]  #"ENSEMBL" 
KEGG_NEED_KEY <- args[7]  #"ENTREZID" 


#scriptPath="/data/script/deseq2+GO+KEGG/Rpipe/"
scriptPath="/Users/Icey/Documents/GitHub/R/"
script1 = paste(scriptPath,"1deseq2-cor-heatmap.R",sep="")
script2 = paste(scriptPath,"2sigGene-volcano.R",sep="")
script3 = paste(scriptPath,"3GO-KEGG.R",sep="")


#1
path1 = "3.DiffExprGene"
dir.create(path1,showWarnings = FALSE)
source(script1)


#2
path2 = "3.DiffExprGene"
dir.create(path2,showWarnings = FALSE)
colData = read.csv(coldata_file, header=T)

for (exp_group in colData$condition[duplicated(colData$condition)])
  if (exp_group != base_group)
    source(script2)

#3
colData = read.csv(coldata_file, header=T)

for (pSet in c(1,0.05)) {
  path3 = paste("4.GO_KEGG",pSet,sep="/")
  dir.create(path3,recursive=TRUE,showWarnings = FALSE)
  for (exp_group in colData$condition[duplicated(colData$condition)]){
    if (exp_group != base_group) source(script3)
  }
}  


# count_table <- "CountMatrix4DESeq.csv"
# coldata_file <- "colData.csv"
# base_group <- "BZ"
# dbname <- 'AH5955'
# kegg_org <- "vda"
# GO_KEY <- "SYMBOL"
# KEGG_NEED_KEY <- "SYMBOL"


