#!/usr/bin/env Rscript

args<-commandArgs(T)

library(dplyr)
library(ggplot2)

plus_file <- args[1]
minus_file <- args[2]
out_file <- args[3]

#read data
df = read.table(file = plus_file, sep = "\t", header=F,stringsAsFactors = F,quote="")
df_minus = read.table(file = minus_file, sep = "\t", header=F,stringsAsFactors = F,quote="")

colnames(df)=c("chr","locus","count")
colnames(df_minus)=c("chr","locus","count")

df$strand = "+"
df_minus$strand = "-"

#log2(sum+1)
df$count <- log2(df$count+1)
df_minus$count <- -log2(df_minus$count+1)

df_merge = rbind(df,df_minus)

#select chr 1-22 X Y 
chr_list <- c(seq(1,22),"X","Y")
df_merge <- df_merge[df_merge$chr %in% chr_list,]

#locus MB
df_merge$locus <- df_merge$locus/1000000

head(df_merge)

df_merge$chr <- factor(df_merge$chr , levels =c(seq(1,22),"X","Y") )

#

png(out_file,height = 2000,width = 2000)
p <- ggplot(df_merge,aes(locus,count))+
    geom_area(aes(fill=strand))+theme_light()+
    theme(axis.text.x = element_text(size = 30), 
          axis.text.y = element_text(size = 10), #y轴刻度大小
            legend.text = element_text(size = 30 ),
         axis.title.x = element_text(size = 30),
         axis.title.y = element_text(size = 30),
         strip.text.y = element_text(size = 30),#分面标签
         legend.title = element_text(size=30)) + #legend 
    xlab("chromosome position (Mb)")+
    ylab("Mean of read density(log2)")+
    scale_y_continuous(breaks = c(-20, 0, 20),labels=c("20","0","20")) #仅显示max和min，去掉负号
                       
p <- p+facet_grid( chr ~ . )  #和wrap差别：分面标签在右侧而不是上方
print(p)
dev.off()