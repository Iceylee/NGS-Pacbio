#!/usr/bin/env Rscript

args<-commandArgs(T)

if (length(args)!=3) {
  print ("3 arguments must be supplied:")
  print ("indir genomeList outdir")
  print ("outData/diamond/ ./genome.list outData/anno/")
  stop(call.=FALSE)
} 

indir <- args[1] 
genomeList <- args[2]
outdir <- args[3] 



#library
library(dplyr)

genome.list = read.table(file=genomeList, header=F,stringsAsFactors = F)


#1.Function Result
result <- genome.list
#View(genome.list)
colnames(result) <- "Gene_ID"

anno_list = c("cog","go","kegg","nr","pfam","sp","trembl")
ANNO_list = c("COG","GO","KEGG","NR","Pfam","SwissProt","TrEMBL")

db_list ="Gene_ID"

for (index in c(rep(1:7))){
  i = anno_list[index]
  I = ANNO_list[index]

  filename = paste(indir,"/",i,".anno.out", sep="")
  outfilename = paste(outdir,"/",i,".annotation.txt",sep="")
  
  #if one anno file not exist then skip this one 
  if (!file.exists(filename)) next

  df = read.table(file = filename, sep = "\t", header=F,stringsAsFactors = F,quote="")
  
  if (i == "go" || i == "kegg"){
    colnames(df) <- c("Gene_ID", "Function")
    gene_func <- df
    #KEGG
    df <- df[!df[,2]=="",]
    df <- df[!df[,1]=="",]

    db_list <- c(db_list,I)
  }
  else{
    colnames(df) <- c("Gene_ID",paste(I,"_ID",sep=""),"Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "Function")
    gene_func<- select_(df,"Gene_ID",paste(I,"_ID",sep=""),"Function")
      #注释为空的行删除
    df <- df[!df[,13]=="",] 

    db_list <- c(db_list,paste(I,"_ID",sep=""),paste(I,"_term",sep=""))  
  }
  #ID为空的行删除
  df <- df[!df[,2]=="",]

  #输出
  write.table(df,file = outfilename,sep = "\t",row.names = F,quote = F)
  result <- left_join(result,gene_func,by = "Gene_ID") 
  
}

result[is.na(result)] <- "-"

colnames(result) <- db_list
write.table(result,file=paste(outdir,"/Annotation_Summary.txt",sep=""),na="",row.names = F,quote = F,sep = "\t")



