#update
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("BiocUpgrade")
# biocLite()

#模式生物
#Zebrafish
# library(org.Dr.eg.db)
#mouse
# source("https://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# library(org.Mm.eg.db)  
#human
#library(org.Hs.eg.db)


library(clusterProfiler)
library(DOSE)

rm=list(ls())

library(clusterProfiler)
library(DOSE)

library(dplyr)

require(AnnotationHub)
hub <- AnnotationHub()
query(hub, "pseudomonas putida")
maize <- hub[['AH56235']]

length(keys(maize))
columns(maize)
keys(maize)[1:100]


gene_list <- read.csv("swissprot-trezID.txt", sep = '\t',header = F,stringsAsFactors=F)
gene_id <- gene_list[,3]

##enrichGO
info_go_BP <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keyType = "ENTREZID", 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
info_go_CC <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keyType = "ENTREZID", 
                       ont = "CC", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)
info_go_MF <- enrichGO(gene = gene_id, 
                       OrgDb = maize, 
                       keyType = "ENTREZID", 
                       ont = "MF", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

#output table
write.table(as.data.frame(info_go_BP@result), file="GO_BP_out.txt",quote=F,row.names = F,sep = "\t")
write.table(as.data.frame(info_go_CC@result), file="GO_CC_out.txt",quote=F,row.names = F,sep = "\t")
write.table(as.data.frame(info_go_MF@result), file="GO_MF_out.txt",quote=F,row.names = F,sep = "\t")

dbDisconnect() 
