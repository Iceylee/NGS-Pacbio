library(clusterProfiler)
library(DOSE)

#ensure organism available
#kegg_code scientific_name common_name
search_kegg_organism('vda', by='kegg_code')

gene_lists <- read.table("Down.Gene.List.txt", sep = '\t',header = F,stringsAsFactors = F)
gene_lists2 <- (gene_lists[,1])

ids <- bitr(gene_lists2, fromType="UNIGENE", toType=c("SYMBOL","ENTREZID"), OrgDb=maize)

gene_id <- ids[,3]
#head(ids)
##NEED SYMBOL 

kk <- enrichKEGG(gene = gene_id, organism = "zma", keyType = "kegg", pvalueCutoff = 0.05)

kk_df <- as.data.frame(kk) %>%
  dplyr::select(-ID)

write.table(kk_df, file="KEGG_out.txt",quote = F,sep = "\t")

##plot
pdf(file="KEGG_dotplot.pdf")
dotplot(kk)
dev.off()


