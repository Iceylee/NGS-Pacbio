setwd("/Users/Icey/work/2018/plot/1012BSA-snpindex/")
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

############
##1.ED plot
############

manhattan_plot <- function(table){
    
    colnames(table)[1:3] <- c("chrom","pos","data")
    
    mycolors = c(brewer.pal(name="Set1", n = 8), brewer.pal(name="Paired", n = 12)) #max 20 colors
    
    img <- ggplot(table, aes(x = pos,y = data,col = chrom)) + 
    geom_point(size=0.3) +
    facet_grid(. ~ chrom,space="free",scale="free",switch="both") + #right label down
    scale_colour_manual(values = mycolors)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), #remove x tick
          legend.position="none",
         strip.background =element_rect(fill="white"),  #remove legend
          panel.background = element_blank(),
          panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line = element_line(), # add axis line
          panel.spacing = unit(0, "lines"),
         plot.title = element_text(hjust = 0.5))+ #remove panel margin 
        geom_smooth(se=F,color="black",method="loess",size=0.5,span=0.1) #span 
    return(img)
}

df <- read.table("Variant_SNP_Result_ED.txt",header=T,sep="\t",stringsAsFactors=F)

chr_list <- c(sprintf("A%02d", 1:10),sprintf("B%02d", 1:10))

ED_data <- select(df,CHROM,POS,ED) %>%
    filter(CHROM %in% chr_list) %>%
    mutate(ED=ED^4)

p <- manhattan_plot(ED_data) + ylab("Euclidean Distance") + xlab(NULL)


#从ggplot提取拟合数据 求阈值 median+3SD
loess_data <- ggplot_build(p)$data[[2]]

threshold_value <- median(loess_data$y) + 3 * sd(loess_data$y)

# 求区间
#所有染色体数据量太大 仅跑A05
ED_A05 <- filter(ED_data,CHROM=="A05")
smooth_vals <- predict(loess(ED~POS,ED_A05), ED_A05$POS) 
result <- cbind(ED_A05,smooth_vals)
confidence_area <- result[result$smooth_vals > threshold_value,] 
range(confidence_area$POS)

#plot
ED_plot <- p + geom_hline(yintercept = threshold_value, lty=2,col="pink",lwd=1)
ggsave("ED.pdf",ED_plot, width=12, height=5, units="in")
ggsave("ED.png",ED_plot, width=12, height=5, units="in")

threshold_value


#2.SNP_index rerun
# df <- read.csv("Finnal.txt",header=T)

# chr_list <- c(sprintf("A%02d", 1:10),sprintf("B%02d", 1:10))
# SNP_index <- filter(df,Chrom %in% chr_list)

# SNP_index1 <- select(SNP_index,Chrom,Pos,FP.5.like_SNP.Index) 
# SNP_index2 <- select(SNP_index,Chrom,Pos,SP.3.like_SNP.Index)
# SNP_index3 <- select(SNP_index,Chrom,Pos,detal_snp_index)

# p1 <- manhattan_plot(SNP_index1)+ ylab("FP_5_like SNP_Index") +xlab(NULL)
# p2 <- manhattan_plot(SNP_index2)+ ylab("SP_3_like SNP_Index") +xlab(NULL)
# p3 <- manhattan_plot(SNP_index3) + ylab("delta SNP_Index") +xlab(NULL)

# p <- arrangeGrob(p1, p2,p3,ncol = 1, heights=c(3,3,3), widths=12)
# ggsave("SNP_index.pdf",p, width=12, height=10, units="in")
# ggsave("SNP_index.png",p, width=12, height=10, units="in")


############
#3. SNP_index & confidence value
############
# install devtools first to download packages from github
#install.packages("devtools")

# use devtools to install QTLseqr
# devtools::install_github("bmansfeld/QTLseqr")

library(QTLseqr)

#data preparation

rawData <- read.csv("Finnal.txt",header=T)

df = rawData
df$CHROM=rawData$Chrom
df$POS=rawData$Pos
df$REF=rawData$Ref
df$ALT=rawData$Alt
df$AD_REF.FP5=rawData$FP.5.like_FP.5_count
df$AD_ALT.FP5=rawData$FP.5.like_SP.3_count
df$AD_REF.SP3=rawData$SP.3.like_FP.5_count
df$AD_ALT.SP3=rawData$SP.3.like_SP.3_count

df <- df[,c(15:22)]
write.table(df,"snp.tsv",sep=",",row.names=FALSE,quote=FALSE)

#begin analysis
df_snp <- importFromTable(file = "snp.tsv", highBulk = "FP5" ,
 lowBulk = "SP3",
chromList = c(sprintf("A%02d", 1:10),sprintf("B%02d", 1:10)))

df_filt <- runQTLseqAnalysis(df_snp, 
                             windowSize = 1e6,
                             popStruc = "F2",
                             depth = 1:100, #need change
                             bulkSize = 20, replications = 10000,
                             intervals = c(95, 99) )

# Ploting
# snp index
df1 <- select(df_filt,CHROM,POS,SNPindex.HIGH)
df2 <- select(df_filt,CHROM,POS,SNPindex.LOW)
p1 <- manhattan_plot(df1)+ ylab("SNP_Index") +xlab(NULL) + labs(title="F bulk")
p2 <- manhattan_plot(df2)+ ylab("SNP_Index") +xlab(NULL) + labs(title="S bulk")
# delta SNP + confidence curve
df3 <- select(df_filt,CHROM,POS,deltaSNP,CI_95,CI_99)
p3 <- manhattan_plot(df3) + ylab("delta SNP_Index") +xlab(NULL) +
    geom_line(aes(y=CI_95),color="green") +geom_line(aes(y=-CI_95),color="green")+
geom_line(aes(y=CI_99),color="red")+geom_line(aes(y=-CI_99),color="red")+labs(title="F-S bulk")
# merge three
p <- arrangeGrob(p1, p2,p3,ncol = 1, heights=c(3,3,3), widths=12)
ggsave("SNP_index.pdf",p, width=12, height=10, units="in")
ggsave("SNP_index.png",p, width=12, height=10, units="in")


# loess拟合：根据POS和deltaSNP做函数，然后根据POS预测deltaSNP 
data_A05 <- filter(data4plot,CHROM=="A05")
smooth_vals <- predict(loess(deltaSNP~POS,data_A05), data_A05$POS) 
result <- cbind(data_A05,smooth_vals)
confidence_area <- result[result$smooth_vals > -(result$CI_99),] 
range(confidence_area$POS)

write.table(df_filt,file="snp_index.tsv",sep="\t",quote=F,row.names=F)

