
library(clusterProfiler)
library(DOSE)
library(dplyr)
require(AnnotationHub)
library(tidyr)

hub <- AnnotationHub()
db <- hub[[dbname]]

groups = paste(exp_group,base_group,sep="vs")
sig_path = paste(path2,"/",groups,"_sig_genes_exprData.txt",sep="")
gene_list <- read.csv(sig_path, sep = '\t',header = T,stringsAsFactors=F)
countData = read.table(count_table, header=TRUE, sep=",")
gene_table = countData[,c(1,2)]
colnames(gene_table) = c("key","id")
gene_list$id = as.character(gene_list$id)
gene_table$id = as.character(gene_table$id)
gene_list2 = left_join(gene_list,gene_table,by="id")
gene_id = gene_list2$key

########################################################
#####################6.GO-result########################
########################################################


##enrichGO
info_go_BP <- enrichGO(gene = gene_id, 
                       OrgDb = db, 
                       keyType = GO_KEY, 
                       ont = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = pSet, 
                       qvalueCutoff = pSet)
info_go_CC <- enrichGO(gene = gene_id, 
                       OrgDb = db, 
                       keyType = GO_KEY, 
                       ont = "CC", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = pSet, 
                       qvalueCutoff = pSet)
info_go_MF <- enrichGO(gene = gene_id, 
                       OrgDb = db, 
                       keyType = GO_KEY, 
                       ont = "MF", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff = pSet, 
                       qvalueCutoff = pSet)

ego_BP <- as.data.frame(info_go_BP@result)
ego_CC <- as.data.frame(info_go_CC@result)
ego_MF <- as.data.frame(info_go_MF@result)



########################################################
#####################7.GO-plot##########################
########################################################

ego_MF_df <- ego_MF %>%
  mutate(onco="Molecular function")
ego_BP_df <- ego_BP %>%
  mutate(onco="Biological process")
ego_CC_df <- ego_CC %>%
  mutate(onco="Cellular component")
ego_three <- rbind(ego_BP_df, ego_CC_df, ego_MF_df)

#top 15 select
ego_three <- ego_three %>%
  arrange_(~ pvalue) %>%
  group_by_(~ onco) %>%
  do(head(., n = 15)) %>%
  arrange(onco,Count)


ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)
ego_three$onco<- factor(ego_three$onco, order=TRUE)
levels(ego_three$onco) <- c("BP","CC","MF")

##plot bar
library(ggplot2)
lable_name <- ego_three$onco[!duplicated(ego_three$onco)]
if (exists("p")) rm(p)
if (dim(ego_three)[1] != 0) {
  p <- ggplot(ego_three, aes(y = Count, x = Description)) +
    geom_bar(stat = "identity", aes(fill = onco), alpha = 1) +
    facet_grid(onco ~ ., scales = "free", space = "free",margins = F) +
    coord_flip()  +
    #scale_y_continuous(limits = c(0, 70))+
    scale_fill_discrete(name = "Ontology", labels = lable_name) +
    theme_light() +
    theme(axis.text = element_text(size = 9), legend.text = element_text(size = 8)) +
    labs(y = "Number of Genes", x = "Term")+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
    scale_y_continuous(labels = function (Count) floor(Count))
  # pdf(file=paste(groups,"_GO_barplot.pdf",sep=""))
  # p
  # dev.off()
  filename  = paste(groups,"_GO_barplot.pdf",sep="")
  ggsave(file=paste(path3,filename,sep="/"),p, width=10, height=10, units="in")
}


########################################################
#####################8.KEGG-plot########################
########################################################
#change ID type
if (GO_KEY != KEGG_NEED_KEY){
  ids <- bitr(gene_id, fromType=GO_KEY, toType=KEGG_NEED_KEY, OrgDb=db)
  id_kegg <- ids[,2]
  filename = paste(groups,"_ID_type.csv",sep="")
  write.csv(ids, paste(path3,filename,sep="/"),row.names = F,quote = F)
  
  }else 
  id_kegg <- gene_id
  
#kegg
kk <- enrichKEGG(gene = id_kegg, organism = kegg_org, keyType = "kegg", pvalueCutoff = pSet)
kk_df <- as.data.frame(kk) %>%
  dplyr::select(-ID)


##plot
if (exists("p")) rm(p)
p = dotplot(kk) + guides(
  color = guide_colorbar(order = 1),
  fill = guide_legend(order = 0))
filename=paste(groups,"_KEGG_dotplot.pdf",sep="")
ggsave(file=paste(path3,filename,sep="/"),p, width=10, height=10, units="in")

###
#output table
ID_conversion <- function(mydf,ID_tab){
  if (nrow(mydf)==0) return(mydf)
  ID_list <- ID_tab[,c(1,2)]
  colnames(ID_list) <- c("geneID","geneName")
  ID_list$geneID=as.character(ID_list$geneID)
  mydf2 <- mydf %>% separate_rows(geneID,sep="/")
  mydf3 <- left_join(mydf2,ID_list,by="geneID")
  mydf4 <- aggregate(mydf3[,c(2:ncol(mydf3))], by = list(mydf3[,1]), function(x) paste(unique(x),collapse="/"))
  mydf5 <- mydf4[,-ncol(mydf4)+2] #delete original ID col
  colnames(mydf5)[1] <- "ID"
  return(mydf5)
}

#GO

ID_name_list = read.table("CountMatrix4DESeq.csv",header=T,sep = ",",stringsAsFactors = F,quote="") 

ego_BP_conv = ID_conversion(ego_BP,ID_name_list)
ego_CC_conv = ID_conversion(ego_CC,ID_name_list)
ego_MF_conv = ID_conversion(ego_MF,ID_name_list)


#output table
fileBP = paste(groups,"_GO_BP_out.txt",sep="")
fileCC = paste(groups,"_GO_CC_out.txt",sep="")
fileMF = paste(groups,"_GO_MF_out.txt",sep="")
write.table(ego_BP_conv, file=paste(path3,fileBP,sep="/"),quote=F,row.names = F,sep = "\t")
write.table(ego_CC_conv, file=paste(path3,fileCC,sep="/"),quote=F,row.names = F,sep = "\t")
write.table(ego_MF_conv, file=paste(path3,fileMF,sep="/"),quote=F,row.names = F,sep = "\t")

#kegg
kegg = kk_df

filename = paste(groups,"_ID_type.csv",sep="")
keggcsvFile=paste(path3,filename,sep="/")


if (file.exists(keggcsvFile)) {
  ID1_ID2_list = read.table(keggcsvFile,header=T,sep = ",",stringsAsFactors = F,quote="") 

  #ID_name_list2=ID_name_list
  ID_name_list2 = read.table("WR180218S_VS_WR180221S.diff.sig",header=T,sep = "\t",stringsAsFactors = F,quote="") 

  colnames(ID_name_list2)[1]="geneID" 
  colnames(ID1_ID2_list)[1]="geneID" 

  kegg_ID_list = left_join(ID1_ID2_list,ID_name_list2,by="geneID")
  kegg_ID_list = kegg_ID_list[,c(2,3)]
} else{
  kegg_ID_list = ID_name_list
}
#if kegg no result do no convert 
# if (nrow(kegg)==0){
#   kegg_conv=kegg
# }else{
# kegg_conv = ID_conversion(kegg,kegg_ID_list)
# }
kegg_conv = ID_conversion(kegg,kegg_ID_list)

filename=paste(groups,"_KEGG_out.txt",sep="")
write.table(kegg_conv, file=paste(path3,filename,sep="/"),quote = F,sep = "\t")


