直接上代码吧：
rm(list=ls())
setwd("C:/Users/lixia/Desktop/XYL/single cell/single cell/ding/")
DEG=read.table("single cell DGE P C.csv",header = T, sep=",",stringsAsFactors = F,row.names = 1)
probeset=rownames(DEG[abs(DEG[,2])>2 & DEG[,5]<0.01,])
length(probeset)
##加载ID转换的包

library(annotate) # lookUp函数是属于annotate这个包的

## **转换时选用何种数据库需要查下！！！**###转换ID
library(biomaRt)  ###useMart函数
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
tmp<-getBM(attributes=c("entrezgene","hgnc_symbol","ensembl_gene_id"), filters = "hgnc_symbol", 
            values =probeset, mart=ensembl)
###不太明白以下缺少了这么多的genne。。。
###获取c中attributes的信息，用filters进行过滤，values是输入的值，mart是指定数据库。。。
head(tmp)
write.table(tmp, file="DGE C VS P.xls", sep="\t",quote=F)
###entrezgene 为 NCBIgene ID，gene的ID号，一般是entrez ID ,ensemble_gene_idw为enemble数据库ID;
##################
EGID <-tmp$entrezgene
length(unique(EGID))
diff_gene_list <- unique(EGID)

### Go analysis BP CC and MF 分析
library(GOstats)
#####BP
paramsBP<- new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation="org.Hs.eg.db",
             
             ontology="BP", pvalueCutoff=0.01, testDirection = "over")
GO.BP = hyperGTest(paramsBP)
BP<-summary(GO.BP)
htmlReport(GO.BP,file="P C_go BP pvaluecutoff 0.01.html")
head(BP)

####CC
paramsCC<- new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation="org.Hs.eg.db",
               ontology="CC", pvalueCutoff=0.01, testDirection = "over")
GO.cc = hyperGTest(paramsCC)
CC<-summary(GO.cc)
htmlReport(GO.cc,file="P C_go CC pvaluecutoff 0.01.html")
head(CC)
###MF
paramsMF<- new("GOHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation="org.Hs.eg.db",
               
               ontology="MF", pvalueCutoff=0.01, testDirection = "over")
GO.MF = hyperGTest(paramsMF)
MF<-summary(GO.MF)
htmlReport(GO.MF,file="P C_go MF pvaluecutoff 0.01.html")
head(MF)

#then do kegg pathway enrichment !
library(org.Hs.eg.db)
hyperG.params = new("KEGGHyperGParams", geneIds=diff_gene_list, universeGeneIds=NULL, annotation="org.Hs.eg.db",
                    categoryName="KEGG", pvalueCutoff=1, testDirection = "over")
KEGG.hyperG.results = hyperGTest(hyperG.params);
htmlReport(KEGG.hyperG.results, file="kegg.enrichment.html", summary.args=list("htmlLinks"=TRUE))