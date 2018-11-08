
library(dplyr)

rm(list=ls())
#C123_1000_gene_MF_GO.txt
# a = "C123_1000_gene_"
# a = "Splice_"
a = "GO"

ego_MF <-read.table(file=paste(a,"_MF_OUT.txt",sep=""),header=T,sep = "\t",stringsAsFactors = F) 
ego_CC<-read.table(file=paste(a,"_CC_OUT.txt",sep=""),header=T,sep = "\t",stringsAsFactors = F)
ego_BP<-read.table(file=paste(a,"_BP_OUT.txt",sep=""),header=T,sep = "\t",stringsAsFactors = F) 


ego_MF_df <- ego_MF %>%
  mutate(onco="Molecular function")
ego_BP_df <- ego_BP %>%
  mutate(onco="Biological process")
ego_CC_df <- ego_CC %>%
  mutate(onco="Cellular component")

ego_three <- rbind(ego_MF_df, ego_CC_df, ego_BP_df)

#top 15 select
#wt = qvalue|Count
ego_three<- ego_three %>%
  group_by(onco) %>%
  top_n(n = 15, wt = Count) %>%
  arrange(onco,Count)

ego_three$Description<- factor(ego_three$Description, order=TRUE, levels=ego_three$Description)
#ego_three$onco<- factor(ego_three$onco, order=TRUE, levels=c("Molecular function","Cellular component","Biological process"))
ego_three$onco<- factor(ego_three$onco, order=TRUE)
levels(ego_three$onco) <- c("BP","CC","MF")

View(ego_three)


##plot bar
library(ggplot2)

#lable_name <- c("Biological process", "Molecular function", "Cellular component")
lable_name <- ego_three$onco[!duplicated(ego_three$onco)]

p <- ggplot(ego_three, aes(y = Count, x = Description)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 1) +
  facet_grid(onco ~ ., scales = "free", space = "free",margins = F) +
  coord_flip()  +
  #scale_y_continuous(limits = c(0, 70))+
  scale_fill_discrete(name = "Ontology", labels = lable_name) +
  theme_light() +
  theme(axis.text = element_text(size = 9), legend.text = element_text(size = 8)) +
  labs(y = "Number of Genes", x = "Term")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) 

# png(file="GO.png", bg="transparent") 
# p
# dev.off()
pdf(file="GO_barplot.pdf")
p
dev.off()

ggsave("GO_barplot.pdf",p, width=10, height=10, units="in")
