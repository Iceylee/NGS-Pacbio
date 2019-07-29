#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=4) {
  print ("4 arguments must be supplied:")
  print ("SigFileDir KEGGDir GeneKeyType Species")
  print ("./nonModel/3.DiffExprGene/ ./nonModel/4.GO_KEGG_Enrichment/0.05/ ENSEMBL ko")
  stop(call.=FALSE)
} 


###Arguments
SIG_PATH = args[1] #"./nonModel/3.DiffExprGene/"
KEGG_DIR = args[2] #"./nonModel/4.GO_KEGG_Enrichment/0.05/"
GENE_KEY= args[3] #"ENSEMBL" 
SPECIES = args[4] #"hsa" or "ko"


suppressMessages(library(pathview))
#############################################
######For Model Animal(eg.Human)#############
#############################################
PathViewModel <- function(sigfile){  
  workdir=getwd()
  print(sigfile)
  #folddir
  groups = sapply(strsplit(sigfile, "_sig"), "[", 1)
  pathway_path = paste(KEGG_DIR,"/",groups,"_Pathway/",sep="")
  
  #check if exists，exists then delete
  if (dir.exists(pathway_path)) {unlink(pathway_path,recursive=TRUE)}
  dir.create(pathway_path,showWarnings = FALSE,recursive=T)
  
  #read in data
  Gene <- read.table(paste(SIG_PATH,sigfile,sep="/"), sep="\t", header = T, row.names = 1)        
  path_info = paste(KEGG_DIR,"/",groups,"_KEGG_Enrichment.tmp.txt",sep="")
  Path <- read.csv(path_info, header = T, row.names = 1, sep="\t")
  
  #Gene expr data
  Gene$geneID <- rownames(Gene)
  
  #去掉括号及括号内字符
  Gene$geneID = gsub("\\s*\\([^\\)]+\\)","",Gene$geneID)
  
  #去除重复
  Gene <- Gene[!duplicated(Gene$geneID),]
  row.names(Gene) <- Gene$geneID
  Gene4Path <- dplyr::select(Gene,log2FoldChange)
  
  #KEGG path id
  KEGG_path_id <- (rownames(Path))
  setwd(pathway_path)

  cc <- function(KEGG_path_id){
    if (num >= 20) {return(NULL)}
    pv.out <- pathview(gene.data = Gene4Path[,1,drop=FALSE], pathway.id = KEGG_path_id, 
                       species = SPECIES, kegg.native = TRUE,
                       gene.idtype = GENE_KEY)
    if (class(pv.out)=="list") {num <<- num + 1}
  }
  
  sapply(KEGG_path_id, cc) 
  setwd(workdir)
  
}

#############################################
######For Non-Model Animal##################
#############################################
PathViewNoModel <- function(sigfile){  
  workdir=getwd()
  #folddir
  groups = sapply(strsplit(sigfile, "_sig"), "[", 1)
  pathway_path = paste(KEGG_DIR,"/",groups,"_Pathway/",sep="")

  #check if exists，exists then delete
  if (dir.exists(pathway_path)) {unlink(pathway_path,recursive=TRUE)}
  dir.create(pathway_path,showWarnings = FALSE,recursive=T)
  
  #read in data
  Gene <- read.table(paste(SIG_PATH,sigfile,sep="/"), sep="\t", header = F,stringsAsFactors=F)        
  path_info = paste(KEGG_DIR,"/",groups,"_KEGG_Enrichment.tmp.txt",sep="")
  Path <- read.csv(path_info, header = T, row.names = 1, sep="\t")
  
  #colnames
  #如果第一行是行名(某列为log2FoldChange)则去掉
  if ("log2FoldChange" %in% Gene[1,]) {
    Gene <- Gene[-1,]
  }
  colnames(Gene)[1:4] = c("koID","geneID","baseMean","log2FoldChange")
  
  #去除重复
  Gene <- Gene[!duplicated(Gene$koID),]
  row.names(Gene) <- Gene$koID
  Gene4Path <- dplyr::select(Gene,log2FoldChange)
  
  #KEGG path id
  KEGG_path_id <- (rownames(Path))
  setwd(pathway_path)
  num = 0
  cc <- function(KEGG_path_id){
    if (num >= 20) {return(NULL)}
    pv.out <- pathview(gene.data = Gene4Path[,1,drop=FALSE], pathway.id = KEGG_path_id, 
                       species = SPECIES, kegg.native = TRUE,
                       gene.idtype = GENE_KEY)
    if (class(pv.out)=="list") {num <<- num + 1}

  }
  
  sapply(KEGG_path_id, cc) 
  setwd(workdir)
  
}


#############################################
######Run Run Run############################
#############################################


if (SPECIES != "ko") {
  sigfiles = list.files(SIG_PATH,pattern="sig_genes_exprData.txt")
  sapply(sigfiles, PathViewModel)
}

if (SPECIES == "ko") {
  sigfiles = list.files(SIG_PATH,pattern="sig_genes_exprData_pathview.txt")
  sapply(sigfiles, PathViewNoModel)
}
