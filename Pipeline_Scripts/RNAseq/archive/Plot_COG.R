#!/usr/bin/env Rscript

args<-commandArgs(T)

if (length(args)!=2) {
  print ("2 arguments must be supplied:")
  print ("countFile plotPrefix")
  print ("COG_Annotation.txt ./cog_bar")
  stop(call.=FALSE)
} 


COGFile <- args[1] #"ED_Genes_COG_Annotation.txt"
plotPrefix <- args[2] #"COG_barplot"

# library
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))

# load data
cog_data <- read_tsv(COGFile,col_names=FALSE)

cog_data <- cog_data[,c(1,4,6)]
colnames(cog_data) <- c("GeneID","class","term")
cog.count <- cog_data %>% group_by(class,term) %>% count()


cog.count$term <- paste(cog.count$class,cog.count$term,sep = " : ")

cog.count$term <- as.factor(cog.count$term)

#term列转为factor后，顺序默认为字母顺序。改为原来的顺序。
#cog.count$term <- factor(cog.count2$term, levels = cog.count2$term)
#因为term和ABC相连接后，字母顺序就是我们要的。因此不用再次排序


#bar plot
p <- ggplot(cog.count,aes(x = class, y = n)) +
  #fill用term列，这样legend就直接打出term了
  geom_bar(aes(fill = term),stat= "identity")+
  #去掉网格和背景颜色。调整lengend标签文件大小
  theme_light() +theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),legend.text=element_text(size=8)) +
  #三个标签 x轴 y轴 legend（fill）
  labs(y = "Number of matched genes", x = "Function class",fill = "Function Classification") +
  #颜色
  scale_color_manual(values=rainbow(25)) +
  #legend从默认的两列调整为一列。 symbol size small
  guides(fill=guide_legend(ncol=1,override.aes = list(size = 0.2)))

#调整legend symbol的大小。 
#   guides(fill=guide_legend(override.aes = list(width = 0.5))) +
#   guides(color=guide_legend(override.aes=list(fill=NA))) +
#   theme(legend.key=element_blank()) 

#存储pdf
ggsave(paste(plotPrefix,".pdf",sep=""),p, width=10, height=7, units="in")

#png
ggsave(paste(plotPrefix,".png",sep=""),p, width=10, height=7, units="in")