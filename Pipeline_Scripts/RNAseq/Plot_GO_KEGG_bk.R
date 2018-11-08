#!/usr/bin/env Rscript

args <- commandArgs(T)

if(length(args) != 1){
	print('Rscript /data1/script/GO_KEGG_Annotation/Plot_GO_KEGG.R "FC2vsFD4"')
	stop(call.=FALSE)
}

library(ggplot2) 
library(dplyr)
library(stringr)
library(clusterProfiler)

#File <- args[1] # FC2vsFD4.GO.Enrichment_out.txt
Name <- args[1] # FC2vsFD4

GO_Name = paste(Name, "_GO_Enrichment.txt", sep="")
GO_pdf = paste(Name, "_GO_barplot.pdf", sep="")

KEGG_Name = paste(Name, "_KEGG_Enrichment.txt", sep="")
KEGG_pdf = paste(Name, "_KEGG_dotplot.pdf", sep="")

GO_Name

GOFile = read.csv(GO_Name, header = T, sep="\t")

GOFile.P <- subset(GOFile, Classification=="BP")
GOFile.C <- subset(GOFile, Classification=="CC")
GOFile.F <- subset(GOFile, Classification=="MF")

GOFile.C.head = head(arrange(as.data.frame(GOFile.C), pvalue), 10)
GOFile.P.head = head(arrange(as.data.frame(GOFile.P), pvalue), 10)
GOFile.F.head = head(arrange(as.data.frame(GOFile.F), pvalue), 10)

GOFile.three <- rbind(GOFile.C.head, GOFile.F.head, GOFile.P.head)

GOFile.three <- GOFile.three[GOFile.three$pvalue<0.1,]

GOFile.three <- GOFile.three %>%
  arrange_(~ pvalue) %>%
  group_by_(~ Classification) %>%
  do(head(., n = 15)) %>%
  arrange(Classification,Count)



GOFile.three$Label<- sapply(as.character(GOFile.three$Description),function(string) {ifelse (nchar(string)>45, paste(substr(string,1,45),"...",sep=""),string)})
#如果label有重复，第二个加上一个句号，以区别
GOFile.three$Label[duplicated(GOFile.three$Label)]<- paste(GOFile.three$Label[duplicated(GOFile.three$Label)],".",sep="")
GOFile.three$Label <- factor(GOFile.three$Label, order=TRUE,levels=GOFile.three$Label)

GOFile.three$Classification<- factor(GOFile.three$Classification, order=TRUE)
levels(GOFile.three$Classification) <- c("BP","CC","MF")




p = ggplot(GOFile.three, aes(x = Label, y = Count)) + 
geom_bar(stat = "identity", aes(fill = Classification), alpha = 1) + 
facet_grid(Classification ~ ., scales = "free", space = "free", margins = F) + coord_flip() +  
theme_light() + 
theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8))  +  
labs(y = "Number of Genes", x = "Term") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


pdf(file=GO_pdf)
p
dev.off()


KEGGFile = read.csv(KEGG_Name, header = T, sep="\t")

tep_num1 = as.character(KEGGFile[2,3])
tep_num2 = strsplit(tep_num1 , "/")
gene_num = as.numeric(tep_num2[[1]][2])

kegg_list <- KEGGFile %>%
  mutate(Ratio = Count/gene_num) %>%
  arrange(pvalue)

kegg_list <- kegg_list[kegg_list$pvalue<0.1,]

kegg_list <- head(kegg_list, n = 10)

kegg_list <- arrange(kegg_list,Ratio)

#View(kegg_list.sorted.head)
p <- ggplot(kegg_list, aes(x = Ratio, y = reorder(Term, Ratio))) + 
  geom_point(aes(size=Count, color=pvalue)) + 
  scale_size("Count") + 
  scale_color_continuous(low="red", high="grey") + 
  theme_light() + 
  labs(x="Gene Ratio", y = "Term") 

pdf(file=KEGG_pdf)
p
dev.off()


