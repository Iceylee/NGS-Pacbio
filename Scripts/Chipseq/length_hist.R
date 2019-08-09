#!/usr/bin/env Rscript

args<-commandArgs(T)

# test if there is at least two argument: if not, return an error
if (length(args)!=3) {
	print ("3 arguments must be supplied:")
	print ("len_table x_max image_prefix")
	print ("1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.seqLengths 3000 length_hist_plot")
  stop(call.=FALSE)
} 

len_table <- args[1] #"Trinity_CD-HIT_0.9.seqLengths" geneID \t length
x_max <- args[2] #"3000"
image_prefix <- args[3] #length_hist_plot

len_df = read.table(len_table,header=F,stringsAsFactors = F,sep="\t")

colnames(len_df)=c("ID","length")

library("ggplot2")
p <- ggplot(data=len_df, aes(len_df$length)) + 
  geom_histogram(breaks = seq(0,x_max,100),
                 fill="blue",
                 col="black")+
  labs(x = "Length(bp)", y = "Count Numbers")

ggsave(paste(image_prefix,".pdf",sep=""),p)

ggsave(paste(image_prefix,".png",sep=""),p)



