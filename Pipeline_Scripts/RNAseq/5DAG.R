#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  print ("4 arguments must be supplied:")
  print ("goAnnoFile sigPath outputPath MODEL")
  print ("./5.Annotation/eggnog.emapper.annotations ./4.DiffExprGene/DESeq2/3.DiffExprGene/ ./6.EnrichmentAnalysis/ TRUE")
  stop(call.=FALSE)
} 


goAnnoFile <- args[1] #eggnog文件路径+文件名
sigPath <- args[2] #差异文件所在路径
outputPath <- args[3] #DAG图的输出目录
MODEL <- args[4] #模式生物否










#source("https://bioconductor.org/biocLite.R")
#chooseBioCmirror()
#biocLite("GO.db")
#biocLite("Rgraphviz")

# load libraries
library(topGO,quietly=T, warn.conflicts=F)
library(GO.db,quietly=T, warn.conflicts=F)
library(Rgraphviz,quietly=T, warn.conflicts=F)
library(stringr,quietly=T, warn.conflicts=F)

options(warn=-1)

dir.create(paste(outputPath,"/DAG_plots/",sep=""),showWarnings = FALSE)

# Function
plot_DAG <- function(sigFile,ONCO){
    
    #group
    groups = sapply(strsplit(sigFile, "_sig"), "[", 1)
    
    # load data
    diffgene <- read.csv(paste(sigPath,"/",sigFile,sep=""), sep="\t", header = T)

    # prepare data
    
    mapFile <- goAnno[,2] 
    names(mapFile) <- gsub(" ", "", goAnno[, 1]) ## 1st col as rownames
    geneID2GO <- lapply(mapFile, function(x) gsub(" ", "", strsplit(x, split = ",")[[1]])) ## split the IDs

    geneNames <- names(geneID2GO)
    #interesting_genes <- str_split_fixed(diffgene$id, "::", 3)[,2]

    geneList <- factor(as.integer(geneNames %in% diffgene[,1]))
    names(geneList) = geneNames

    ## create GOData object
                        
    GOdata <- new("topGOdata", ontology=ONCO, allGenes = geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize=5)

    ## run analysis
                        
    resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    allRes <- GenTable(GOdata, classic = resultFis, orderBy = "weight", ranksOf = "classic", topNodes = 20)
    #showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
                        
    workdir=getwd()
    setwd(paste(outputPath,"/DAG_plots/",sep=""))
    printGraph(GOdata, resultFis, firstSigNodes = 5, fn.prefix = paste(groups,"_",ONCO,sep=""), useInfo = "all", pdfSW = TRUE)
    setwd(workdir)
                    
}

#################
###RUN RUN RUN###
#################

sigfiles = list.files(sigPath,pattern="sig_genes_exprData.txt")

                        
if (MODEL == TRUE) {
    system(paste("cat ",goAnnoFile,"|cut -f 1,2|awk -f /data1/script/deseq2+GO+KEGG/Rpipe/merge.awk >",outputPath,"/Gene2GO.tsv.temp",sep=""))
    goAnno <- read.table(paste(outputPath,"/Gene2GO.tsv.temp",sep=""),sep="\t",comment="#",header=F,quote="",stringsAsFactors=F)
}
                        
if (MODEL == FALSE) {
    
    eggnog <- read.table(goAnnoFile,sep="\t",comment="#",header=F,quote="",stringsAsFactors=F)
    goAnno <- eggnog[,c(1,6)] ##1st Gene ID;6th GO col;no col names
}                        
                        
                        
for (i in c("BP","CC","MF")){
    sapply(sigfiles,plot_DAG,ONCO=i)
}
