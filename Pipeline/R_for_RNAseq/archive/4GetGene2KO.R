#!/usr/bin/env Rscript

# This script will convert gene symbols to kegg (KO) identifiers

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=3) {
  print ("3 arguments must be supplied:")
  print ("dbName keggName outpath")
  print ("AH134 hsa ./4.GO_KEGG_Enrichment/")
  stop(call.=FALSE)
} 

dbname <- args[1] #"AH59553"# 'AH57973'
kegg_org <- args[2]  #"vda" # "hsa"
output_path <- args[3] #"./"

suppressMessages(library(DOSE))
suppressMessages(library(dplyr))
suppressMessages(library(AnnotationHub))
suppressMessages(library(clusterProfiler))
suppressMessages(library(RCurl))
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))

hub <- AnnotationHub()
db <- hub[[dbname]]

#提取物种所有Symbol ID
gene_list <- keys(db,keytype="SYMBOL")

#del
gene_list = gene_list[1:10]
#del

# For each gene, look up the identifiers
KEGG_ID = c()
for (g in 1:length(gene_list)) {
  query = paste("http://rest.kegg.jp/find/",kegg_org,"/",gene_list[g],sep="")
  result = getURL(query)
  KEGG_ID = c(KEGG_ID,strsplit(result,"\t")[[1]][1])  
}

# SYMBOL 与 KEGG_ID(hsa:1111) 对应
gene_list = cbind(gene_list,KEGG_ID)
gene_list = as.data.frame(gene_list)

#去掉没有结果的 
gene_list = gene_list[grep(kegg_org,KEGG_ID),]
colnames(gene_list) = c("SYMBOL","KEGG_ID")

# Now look up KO ID
ko = c()
for (g in 1:length(gene_list$KEGG_ID)) {
  query = paste("http://rest.kegg.jp/link/ko/",as.character(gene_list$KEGG_ID[g]),sep="")
  result = getURL(query) 
  
  ko_ko_number = as.character(gsub("\n","",strsplit(result,"\t")[[1]][2]))
    
  if (!is.na(ko_ko_number)){
  ko_number = strsplit(ko_ko_number,":")[[1]][[2]]
  }else {ko_number = ko_ko_number}
  ko = c(ko,ko_number)  
}

# combine
ko = as.character(ko)
df = cbind(gene_list,ko)

ids <- bitr(df$SYMBOL, fromType="SYMBOL", toType="ENTREZID", OrgDb=db)

df <- left_join(df,ids,by="SYMBOL")

df <- df[,c("SYMBOL","ENTREZID","KEGG_ID","ko")]

# Write to file
outfile = paste(output_path,"AllGene_KEGG_Annotation.txt",sep="")
write.table(df,file=outfile,sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)


