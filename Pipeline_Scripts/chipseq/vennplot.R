#!/usr/bin/env Rscript

args<-commandArgs(T)

smp1_stats=as.numeric(args[1])
smp2_stats=as.numeric(args[2])
comm_stats=as.numeric(args[3])
smp1=args[4]
smp2=args[5]

library(VennDiagram)

png(file=paste(smp1,"vs",smp2,"_vennplot.png",sep=""))
draw.pairwise.venn(smp1_stats+comm_stats, smp2_stats+comm_stats, comm_stats, 
                   category = c(smp1, smp2), lty = rep("blank", 
    2), fill = c('#0099CC','#FF6666'), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.025, 2))
dev.off()