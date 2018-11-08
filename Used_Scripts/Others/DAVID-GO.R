ego_CC <- read.csv("result_out/GO-CC.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ego_MF <- read.csv("result_out/GO-MF.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ego_BP <- read.csv("result_out/GO-BP.txt", sep = "\t", header = TRUE,stringsAsFactors = FALSE)
View(ego_CC)
View(ego_MF)
View(ego_BP)

#term 转换 
temp <- sapply(ego_CC$Term, function(x) (strsplit(x,"~"))[[1]])
ego_CC$Term <- temp[2,]

temp <- sapply(ego_BP$Term, function(x) (strsplit(x,"~"))[[1]])
ego_BP$Term <- temp[2,]
temp <- sapply(ego_MF$Term, function(x) (strsplit(x,"~"))[[1]])
ego_MF$Term <- temp[2,]

#p value top 15 select
ego_CC_df <- ego_CC[1:7, c(1:3, 7)] %>%
  transform(onco = "Cellular component", GeneRatio = Count / List.Total * 100)

ego_BP_df <- ego_BP[1:7, c(1:3, 7)] %>%
  transform(onco = "Biological process", GeneRatio = Count / List.Total * 100)

ego_MF_df <- ego_MF[1:7, c(1:3, 7)] %>%
  transform(onco = "Molecular function", GeneRatio = Count / List.Total * 100)

# View(ego_MF_df)
# View(ego_BP_df)
# View(ego_CC_df)

##bind three onco
ego_three <- rbind(ego_CC_df, ego_BP_df, ego_MF_df)
# View(ego_three)

##plot bar
#########NEED CHANGE limits,aes(y,x)
library(ggplot2)

p <- ggplot(ego_three, aes(y = GeneRatio, x = Term)) +
  geom_bar(stat = "identity", aes(fill = onco), alpha = 0.7) +
  facet_grid(onco ~ ., scales = "free") +
  coord_flip()  +
  scale_y_continuous(limits = c(0, 8))+
  scale_fill_discrete(name = "Ontology", labels = c("Biological process", "Molecular function", "Cellular component")) +
  theme_light() +
  theme(axis.text = element_text(size = 8, face = "bold"), legend.text = element_text(size = 8, face = "bold")) +
  labs(y = "Percentage of Genes", x = "Term")

pdf(file="GO_barplot.pdf")
p
dev.off()
