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

GO_Name = paste(Name, ".GO.Enrichment_out.txt", sep="")
GO_pdf = paste(Name, "_GO_barplot.pdf", sep="")

KEGG_Name = paste(Name, ".KEGG.Enrichment_out.txt", sep="")
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

GOFile.three <- GOFile.three %>%
  arrange_(~ pvalue) %>%
  group_by_(~ Classification) %>%
  do(head(., n = 15)) %>%
  arrange(Classification,Input.number)



GOFile.three$Description<- sapply(as.character(GOFile.three$Description),function(string) {ifelse (nchar(string)>45, paste(substr(string,1,45),"...",sep=""),string)})

GOFile.three$Description<- factor(GOFile.three$Description, order=TRUE, levels=GOFile.three$Description)
GOFile.three$Classification<- factor(GOFile.three$Classification, order=TRUE)
levels(GOFile.three$Classification) <- c("BP","CC","MF")




p = ggplot(GOFile.three, aes(x = Description, y = Input.number)) + 
geom_bar(stat = "identity", aes(fill = Classification), alpha = 1) + 
facet_grid(Classification ~ ., scales = "free",margins = F) + coord_flip() +  
theme_light() + 
theme(axis.text = element_text(size = 8), legend.text = element_text(size = 8))  +  
labs(y = "Number of Genes", x = "Term") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())


pdf(file=GO_pdf)
p
dev.off()






KEGGFile = read.csv(KEGG_Name, header = T, sep="\t")
kegg_list <- KEGGFile %>%
  mutate(Ratio = Input.number/Background.number) %>%
  arrange(Ratio)

kegg_list <- head(kegg_list, n = 10)

#View(kegg_list.sorted.head)
p <- ggplot(kegg_list, aes(x = Ratio, y = reorder(Term, Ratio))) + 
  geom_point(aes(size=Input.number, color=pvalue)) + 
  scale_size("Count") + 
  scale_color_continuous(low="red", high="grey") + 
  theme_light() + 
  labs(x="Gene Ratio", y = "Term") 

pdf(file=KEGG_pdf)
p
dev.off()

