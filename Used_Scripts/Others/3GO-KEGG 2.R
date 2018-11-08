
library(clusterProfiler)
library(DOSE)
library(dplyr)
require(AnnotationHub)
library(tidyr)

hub <- AnnotationHub()
db <- hub[[dbname]]

groups = paste(exp_group,base_group,sep="vs")
sig_path = paste("result/2sigGene-volcano/",groups,"_sig_genes_exprData.txt",sep="")
gene_list <- read.csv(sig_path, sep = '\t',header = T,stringsAsFactors=F)

countData = read.table(count_table, header=TRUE, sep=",")
gene_table = countData[,c(1,2)]
colnames(gene_table) = c("key","id")
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

#output table
fileBP = paste(groups,"_GO_BP_out.txt",sep="")
fileCC = paste(groups,"_GO_CC_out.txt",sep="")
fileMF = paste(groups,"_GO_MF_out.txt",sep="")
write.table(as.data.frame(info_go_BP@result), file=paste(path3,fileBP,sep="/"),quote=F,row.names = F,sep = "\t")
write.table(as.data.frame(info_go_CC@result), file=paste(path3,fileCC,sep="/"),quote=F,row.names = F,sep = "\t")
write.table(as.data.frame(info_go_MF@result), file=paste(path3,fileMF,sep="/"),quote=F,row.names = F,sep = "\t")





########################################################
#####################7.GO-plot##########################
########################################################

ego_BP<-read.table(file=paste(path3,fileBP,sep="/"),header=T,sep = "\t",stringsAsFactors = F,quote="") 
ego_CC<-read.table(file=paste(path3,fileCC,sep="/"),header=T,sep = "\t",stringsAsFactors = F,quote="")
ego_MF <-read.table(file=paste(path3,fileMF,sep="/"),header=T,sep = "\t",stringsAsFactors = F,quote="") 

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
rm(p)
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
  filename = paste(groups,"_GO_plot_data.txt",sep="")
  write.table(ego_three, file=paste(path3,filename,sep="/"),quote=F,row.names = F,sep = "\t")
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
filename=paste(groups,"_KEGG_out.txt",sep="")
write.table(kk_df, file=paste(path3,filename,sep="/"),quote = F,sep = "\t")


##plot
rm(p)
p = dotplot(kk) + guides(
  color = guide_colorbar(order = 1),
  fill = guide_legend(order = 0))
filename=paste(groups,"_KEGG_dotplot.pdf",sep="")
ggsave(file=paste(path3,filename,sep="/"),p, width=10, height=10, units="in")



