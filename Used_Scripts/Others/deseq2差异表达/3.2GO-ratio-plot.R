
#得到GO结果：info_go_CC info_go_BP info_go_MF

#get gene sum number
#可能的问题：BP结果为空，无法计算
tep_num1 = as.data.frame(info_go_BP)[1,3]
tep_num2 = strsplit(tep_num1 , "/")
gene_num = as.numeric(tep_num2[[1]][2])




#加上onco列并计算gene-ratio
ego_CC_df <- as.data.frame(info_go_CC)%>%
  transform(onco = "Cellular component", GeneRatio = Count / gene_num * 100)

ego_BP_df <- as.data.frame(info_go_BP)%>%
  transform(onco = "Biological process", GeneRatio = Count / gene_num * 100)

ego_MF_df <- as.data.frame(info_go_MF)%>%
  transform(onco = "Molecular function", GeneRatio = Count / gene_num * 100)

#如果某结果为空，则无法bind 
# ego_three <- rbind(ego_BP_df,ego_MF_df)
ego_three <- rbind(ego_CC_df, ego_BP_df, ego_MF_df)

#取pvalue前10
ego_three<- ego_three %>%
  group_by(onco) %>%
  top_n(n = 15, wt = qvalue) %>%
  arrange(onco,GeneRatio)


##plot bar

#lable_name <- c("Biological process", "Molecular function", "Cellular component")
lable_name <- ego_three$onco[!duplicated(ego_three$onco)]

library(ggplot2)
p <- ggplot(ego_three, aes(y = GeneRatio, x = Description)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 0.7) +
  facet_grid(onco ~ ., scales = "free") +
  coord_flip()  +
  scale_fill_discrete(name = "Ontology", labels = lable_name) +
  theme_light() +
  theme(axis.text = element_text(size = 6, face = "bold"), legend.text = element_text(size = 6, face = "bold")) +
  labs(y = "Percentage of Genes", x = "Term")

pdf(file="GO_barplot.pdf")
p
dev.off()

dbDisconnect() 
