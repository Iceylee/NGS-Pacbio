#!/bin/bash



if [  -f $1/$2 ]; then
	file=$1/$2
	outputFile=4.GO_KEGG_Enrichment/AllGene_KEGG_Annotation.txt


	
fi


  #按照tab，将第二列转换回第一列
  #输出geneID tab koPath
  system(sprintf("bash /data1/script/deseq2+GO+KEGG/Rpipe/gene_ko.sh %s %s" ,path3,filename),intern=TRUE)


cat NCvsSH1_KEGG_Enrichment.txt|awk '
BEGIN{FS=OFS="\t"}
{split($8,x,"/");
for (i in x){print x[i],$1}}'|sed '1 d'|sort -k1,1|sed '1 i geneID\tKEGG_ID' > NCvsSH1_gene_ko.txt