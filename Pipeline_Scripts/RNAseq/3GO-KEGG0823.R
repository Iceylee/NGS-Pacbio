#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least one argument: if not, return an error
if (length(args)!=5) {
  print ("5 arguments must be supplied:")
  print ("dbName keggName keyType keggType")
  print ("AH134 hsa ENSEMBL ENTREZID")
  stop(call.=FALSE)
} 

dbname <- args[1] #'AH134'
kegg_org <- args[2]  #"hsa"
GO_KEY <- args[3]  #"ENSEMBL" 
KEGG_NEED_KEY <- args[4]  #"ENTREZID" 
output_path <- args[5] #"./"


suppressMessages(library(clusterProfiler))
suppressMessages(library(DOSE))
suppressMessages(library(dplyr))
suppressMessages(library(AnnotationHub))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))

GOKEGG <- function(file,pSet)  {
  path3 = paste(output_path,"4.GO_KEGG_Enrichment/",pSet,"/",sep="")
  dir.create(path3,showWarnings = FALSE,recursive=T)
  
  sig_path = paste(path1,"/",file,sep="")
  gene_list <- read.table(sig_path,row.names=1,sep = '\t',header = T,stringsAsFactors=F)
  gene_list$geneID <- rownames(gene_list)

  # if first col type: geneID(geneSymbol)
  # extract symbol in parentesis
  if (grepl("\\(",gene_list$geneID[1])){
    gene_list$geneID <- sapply(gene_list$geneID,function(string){regmatches(string, gregexpr("(?<=\\().*?(?=\\))", string, perl=T))[[1]]})
  }

  groups = sapply(strsplit(file, "_sig"), "[", 1)
    
  gene_id = as.character(gene_list$geneID)

  print(gene_id[1])
  
  ########################################################
  #####################6.GO-result########################
  ########################################################
  
  
  ##enrichGO
  info_go_BP <- enrichGO(gene = gene_id, 
                         OrgDb = db, 
                         keyType = GO_KEY, 
                         ont = "BP", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 1,
                         minGSSize = 1)
  info_go_CC <- enrichGO(gene = gene_id, 
                         OrgDb = db, 
                         keyType = GO_KEY, 
                         ont = "CC", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 1,
                         minGSSize = 1)
  info_go_MF <- enrichGO(gene = gene_id, 
                         OrgDb = db, 
                         keyType = GO_KEY, 
                         ont = "MF", 
                         pAdjustMethod = "BH", 
                         pvalueCutoff = 1, 
                         qvalueCutoff = 1,
                         minGSSize = 1)

  #filter pvalue by pSet
  ego_BP <- as.data.frame(info_go_BP@result) %>% filter(pvalue <= pSet)
  ego_CC <- as.data.frame(info_go_CC@result) %>% filter(pvalue <= pSet)
  ego_MF <- as.data.frame(info_go_MF@result) %>% filter(pvalue <= pSet)
  
  fileBP = paste(groups,"_GO_BP_out.txt",sep="")
  fileCC = paste(groups,"_GO_CC_out.txt",sep="")
  fileMF = paste(groups,"_GO_MF_out.txt",sep="")

  write.table(ego_BP, file=paste(path3,fileBP,sep="/"),quote=F,row.names = F,sep = "\t")
  write.table(ego_CC, file=paste(path3,fileCC,sep="/"),quote=F,row.names = F,sep = "\t")
  write.table(ego_MF, file=paste(path3,fileMF,sep="/"),quote=F,row.names = F,sep = "\t")
  
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
  
  ego_three$Description<- sapply(as.character(ego_three$Description),function(string) {ifelse (nchar(string)>45, paste(substr(string,1,45),"...",sep=""),string)})
  
  ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)
  ego_three$onco<- factor(ego_three$onco, order=TRUE)
  levels(ego_three$onco) <- c("BP","CC","MF")
  
  ##plot bar
  lable_name <- ego_three$onco[!duplicated(ego_three$onco)]
  if (exists("p")) rm(p)
  
  if (dim(ego_three)[1] != 0) {
    p <- ggplot(ego_three, aes(y = Count, x = Description)) +
      geom_bar(stat = "identity", aes(fill = onco), alpha = 1) +
      facet_grid(onco ~ ., scales = "free", space = "free",margins = F) +
      coord_flip()  +
      scale_fill_discrete(name = "Ontology", labels = lable_name) +
      theme_light() +
      theme(axis.text = element_text(size = 9), legend.text = element_text(size = 8)) +
      labs(y = "Number of Genes", x = "Term")+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
      scale_y_continuous(labels = function (Count) floor(Count))

    filename  = paste(groups,"_GO_barplot.pdf",sep="")
    ggsave(file=paste(path3,filename,sep="/"),p, width=10, height=10, units="in")

    print("------------------------------")
    print(pSet)
    print("GO_plot OK!")
    print("------------------------------")

  }
  

  
  ########################################################
  #####################8.KEGG-plot########################
  ########################################################
  #change ID type
  if (GO_KEY != KEGG_NEED_KEY){
    ids <- bitr(gene_id, fromType=GO_KEY, toType=KEGG_NEED_KEY, OrgDb=db)
    id_kegg <- ids[,2]
    filename_tab = paste(groups,"_ID_type.tab",sep="")
    write.table(ids, paste(path3,filename_tab,sep="/"),row.names = F,quote = F,sep="\t")
    
    }else 
    id_kegg <- gene_id
  
  #kegg
  kk <- enrichKEGG(gene = id_kegg, organism = kegg_org, keyType = "kegg", pvalueCutoff = 1,qvalueCutoff = 1,minGSSize = 2)
  kk_df <- as.data.frame(kk) %>%
    #dplyr::select(-ID) %>%
    filter(pvalue <= pSet) %>%
    arrange(pvalue)%>%
    do(head(., n = 15)) %>%
    arrange(GeneRatio)

   #gene sum 
  tep_num1 = kk_df$GeneRatio[1]
  tep_num2 = strsplit(tep_num1 , "/")
  gene_num = as.numeric(tep_num2[[1]][2])

  kk_df <- kk_df %>% transform(GeneRatio = Count / gene_num )

  filename=paste(path3,groups,"_KEGG_Enrichment.txt",sep="")
  write.table(kk_df, file=filename,quote = F,sep = "\t",row.names = FALSE)
  
  ##plot
  if (exists("p")) rm(p)
  if (dim(kk)[1] != 0) {
    pp <- ggplot(kk_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) + geom_point(aes(size=Count, color=pvalue)) + 
      scale_size("Count") + 
      scale_color_continuous(low="red", high="grey") + 
      theme_light() + 
      labs(x="Gene Ratio", y = "Term") 
    filename=paste(groups,"_KEGG_dotplot.pdf",sep="")
    ggsave(file=paste(path3,filename,sep="/"),pp, width=10, height=10, units="in")
    
    print("------------------------------")
    print(pSet)
    print("KEGG_plot OK!")
    print("------------------------------")
  }


  if (exists("filename_tab")) system(sprintf("bash /data1/script/deseq2+GO+KEGG/Rpipe/ID_conv0.sh %s %s" ,path3,filename_tab),intern=TRUE)
}



hub <- AnnotationHub()
db <- hub[[dbname]]



path1 = paste(output_path,"3.DiffExprGene",sep="")
sigfiles = list.files(path1,pattern="sig_genes_exprData.txt")

sapply(sigfiles,GOKEGG,pSet=1)

print("------------------------------")
print("########pSet 1 OK!##########")
print("------------------------------")

sapply(sigfiles,GOKEGG,pSet=0.05)

print("------------------------------")
print("########pSet 0.05 OK!#########")
print("------------------------------")



