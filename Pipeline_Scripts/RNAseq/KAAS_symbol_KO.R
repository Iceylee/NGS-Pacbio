## 1.KEGG
#某物种的所有基因（entrezID），用clusterprofiler做KEGG注释，然后针对每个path拆分，输出每个geneID对应的path。

suppressMessages(library(clusterProfiler))
suppressMessages(library(DOSE))
suppressMessages(library(dplyr))
suppressMessages(library(AnnotationHub))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

dbname <-  "AH59553"# 'AH57973'
kegg_org <- "vda" # "hsa"
GO_KEY <- "SYMBOL" 
KEGG_NEED_KEY <- "ENTREZID" 
output_path <- "./"

hub <- AnnotationHub()
db <- hub[[dbname]]



########################################################
#####################KEGG&Gene Mapping #################
########################################################
Mapping <- function()  {
  path3 = paste(output_path,"4.GO_KEGG_Enrichment/",sep="")
  filename=paste(path3,"AllGene_KEGG_Annotation.txt",sep="")

  if (!exists(filename)) {
      
    #提取所有基因
    gene_list <- keys(db,keytype=KEGG_NEED_KEY)
    kk <- enrichKEGG(gene = gene_list, organism = kegg_org, keyType = "kegg", pvalueCutoff = 1,
                   qvalueCutoff = 1,minGSSize = 1,maxGSSize = 100000000)
    #拆分geneID，生成df。每行输出一个geneID和对应的koID
    GeneSplit <- function(geneID,koID){
      gene_list <- unlist(strsplit(geneID,split="/"))
      return (data.frame(geneID=gene_list,ko_ID=rep(koID,length(gene_list))))
    }

    #将每行产生的df合并。
    df <- data.frame(geneID=NULL,ko_ID=NULL)
    for (i in c(1:nrow(kk))){
      df <- rbind(df,GeneSplit(kk$geneID[i],kk$ID[i]))}

    #按geneID排序
    df <- df %>% arrange(geneID)

    #输出
    
    write.table(df, file=filename,quote = F,sep = "\t",row.names = FALSE)
  }


}

Mapping()

## 2.KEGG API
#提取某物种的所有symbol号，通过KAAS得到对应的kegg 识别号，然后对应的相应的大KO号（pathview 非模式生物需要）。
#KEGG识别号和大KO号是一一对应的。

#### 2.1 测试

#通过symbol得到kegg识别号
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))
  query = paste("http://rest.kegg.jp/find/ath/","CYCB3",sep="")
  result = getURL(query)
  result


#通过 kegg 识别号 得到 大K号
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))
  query = paste("http://rest.kegg.jp/link/ko/","hsa:54776",sep="")
  result = getURL(query)
#   ko = c(ko,as.character(gsub("\n","",strsplit(result,"\t")[[1]][2])))  
 result

search_kegg_organism('Arabidopsis thaliana', by='scientific_name')

#### 2.2 脚本

dbname <- 'AH57973'
kegg_org <- "hsa"
output_path <- "./"

suppressMessages(library(DOSE))
suppressMessages(library(dplyr))
suppressMessages(library(AnnotationHub))
suppressMessages(library(clusterProfiler))

hub <- AnnotationHub()
db <- hub[[dbname]]

# This script will convert a list of gene symbols to kegg (KO) identifiers

# We need libraries to read and parse URL
library(RCurl)
options(RCurlOptions = list(cainfo = system.file("CurlSSL", "cacert.pem", package = "RCurl")))

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

####GOOD CODE
# Write to file
outfile = gsub(".csv","_ko.csv",infile)
write.table(data,file=outfile,row.names=FALSE,col.names=TRUE,sep="\t")