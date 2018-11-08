#两个表，共同有swissID
#小表 通过swiss ID 得到 大库的所有其他信息（若id有重复信息，会多列输出）

setwd("~/work/11月/1116go-id/")
library(dplyr)


genome.list = read.table(file="genomeNum_swissprot.id", sep =" ",header=F,stringsAsFactors = F)

go_db = read.table(file="goa_1000.gaf", sep ="\t",header=F,stringsAsFactors = F)


colnames(genome.list)[2] <- "gene_id"
colnames(go_db)[2] <- "gene_id"

go_db_small <- go_db[,c(2,5,9)]

#######
result <- left_join(genome.list,go_db_small,by = "gene_id")


# count import items
# sum(!is.na(result$phi_Function))

write.csv(result,file="result.go.csv",na="",row.names = F)