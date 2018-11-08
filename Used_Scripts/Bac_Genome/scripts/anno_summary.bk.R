setwd("/Users/Icey/work/2018/0105_852Brevibacillus唐兵/original/3功能分析/result/")
library(dplyr)

dir.create("out/")
genome.list = read.table(file="genome.list", header=F,stringsAsFactors = F)


#1.Function Result
result <- genome.list
#View(genome.list)
colnames(result) <- "Gene_ID"

anno_list = c("cog","go","kegg","nr","pfam","sp","trembl")
ANNO_list = c("COG","GO","KEGG","NR","Pfam","SwissProt","TrEMBL")

for (index in c(rep(1:7))){
  i = anno_list[index]
  I = ANNO_list[index]

  filename = paste(i,".anno.out", sep="")
  outfilename = paste("out/",i,".annotation.txt",sep="")
  df = read.table(file = filename, sep = "\t", header=F,stringsAsFactors = F,quote="")
  
  if (i == "go" || i == "kegg"){
    colnames(df) <- c("Gene_ID", "Function")
    gene_func <- df
    #KEGG
    df <- df[!df[,2]=="",]
    df <- df[!df[,1]=="",]
  }
  else{
    colnames(df) <- c("Gene_ID",paste(I,"_ID",sep=""),"Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "Function")
    gene_func<- select(df,Gene_ID,Function)
      #注释为空的行删除
    df <- df[!df[,13]=="",]   
  }
  #ID为空的行删除
  df <- df[!df[,2]=="",]

  #输出
  write.table(df,file = outfilename,sep = "\t",row.names = F,quote = F)
  result <- left_join(result,gene_func,by = "Gene_ID") 
}

result[is.na(result)] <- "-"

colnames(result) <- c("Gene_ID","COG","GO","KEGG","NR","Pfam","SwissProt","TrEMBL")
write.table(result,file="out/Annotation_Summary.txt",na="",row.names = F,quote = F,sep = "\t")


#2.Pathogen Result
anno_list = c("ardb","phi","vfdb","cazy")
ANNO_list = c("ARDB","PHI","VFDB","CAZy")

for (index in c(rep(1:4))){
  i = anno_list[index]
  I = ANNO_list[index]

  filename = paste(i,".anno.out", sep="")
  outfilename = paste("out/",i,".annotation.txt",sep="")
  df = read.table(file = filename, sep = "\t", header=F,stringsAsFactors = F,quote="")
  colnames(df) <- c("Gene_ID",paste(I,"_ID",sep=""),"Identity","Align_Length","Mismatch","Gap","Q_start","Q_end", "T_start", "T_end", "E_value", "Score", "Function")
  #注释为空的行删除
  df <- df[!df[,13]=="",]   
  #ID为空的行删除
  df <- df[!df[,2]=="",]

  #输出
  write.table(df,file = outfilename,sep = "\t",row.names = F,quote = F)
}
