################################################
##expression
################################################
CountTable = /data2/ClientData/2018_06/YangFengLian/project2/R_input/CountMatrix4DESeq1.csv
ColDataFile = /data2/ClientData/2018_06/YangFengLian/project2/R_input/colData.csv
BaseGroup = MD

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R /data2/ClientData/2018_06/YangFengLian/project2/R_input/CountMatrix4DESeq1.csv /data2/ClientData/2018_06/YangFengLian/project2/R_input/colData.csv MB ./

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R /data2/ClientData/2018_06/YangFengLian/project2/R_input/CountMatrix4DESeq1.csv /data2/ClientData/2018_06/YangFengLian/project2/R_input/colData.csv MC ./

cat MDvsMB.GO.Enrichment_out.txt|awk -F"\t" '$9=="BP" {print $1,$2,$3,$5}'
#sort 分隔符 逆序
cat MDvsMB.GO.Enrichment_out.txt|awk 'BEGIN{FS=OFS="\t"} $9=="BP" {print $1,$2,$3,$5}'|sort -t $'\t' -k4n,4|head




################################################
##enrichment
################################################
group=MDvsMC 
#MDvsMB

#去掉 .p1
cat ${group}_sig_genes_exprData.txt|sed "s/\.p.//g" > ${group}_sig_genes_exprData.txt.temp
cd 4.GO_KEGG/

#### KEGG
# 1. Annotation
awk -F'\t' 'NR==FNR{a[$1]=$0;next}{if(a[$1]) {print $0}}' ../3.DiffExprGene/${group}_sig_genes_exprData.txt.temp  ../../Genome/KEGG.annotation.txt > ${group}.KEGG.annotation.txt
# 2. Enrich
python /data1/script/GO_KEGG_Annotation/Pyher_KEGGEnrichment_Analysis.py  ${group}.KEGG.annotation.txt ../../Genome/KEGG.annotation.txt ${group}

#### GO
# 1. Annotation
awk -F'\t' 'NR==FNR{a[$1]=$0;next}{if(a[$1]) {print $0}}' ../3.DiffExprGene/${group}_sig_genes_exprData.txt.temp  ../../Genome/GO.annotation.txt  > ${group}.GO.annotation.txt
# 2. Enrich
python  /data1/script/GO_KEGG_Annotation/Pyher_GOEnrichment_Analysis.py ${group}.GO.annotation.txt ../../Genome/GO.annotation.txt  ${group}

#### Plot
Rscript /data1/script/GO_KEGG_Annotation/Plot_GO_KEGG.R  








