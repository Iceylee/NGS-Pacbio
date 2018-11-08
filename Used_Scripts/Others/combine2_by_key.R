
#Rscript combine2_by_key.R genomeNum_swissprot.id GO.out gene_swiss_GO.id 2 1
#输入文件1 输入文件2 输出文件 key1 key2

Args <- commandArgs()

library(dplyr)


#test
# Args<-vector(mode="character",length=10)
# Args[6:10] <-c("genomeNum_swissprot.id","GO.out","gene_swiss_GO.id","2","1")


input_file1 = Args[6]
input_file2 = Args[7]
output_file = Args[8]
file1_col_num = as.numeric(Args[9]) #2
file2_col_num = as.numeric(Args[10]) #1


table1 = read.csv(file=input_file1, sep ="\t",header=F,stringsAsFactors = F)

table2 = read.csv(file=input_file2, sep ="\t",header=F,stringsAsFactors = F)

colnames(table1)[file1_col_num] <- "key"
colnames(table2)[file2_col_num] <- "key"


#######
result <- dplyr::left_join(table1,table2,by = "key")
result <- result[!duplicated(result$key), ] #去重复 

write.table(result,file=output_file,na="",row.names = F,quote = F,sep = "\t", col.names = F)
