#JiangYang

cat ../Gene2Symbol.txt|awk 'BEGIN{FS="\t";OFS=","}{$1=$1}{print $0}' > Gene2Symbol.csv

awk '
BEGIN{FS=OFS=","}
NR==FNR{a[$1]=$2;next}
NR!=FNR && FNR==1 {print $0}
NR!=FNR && FNR!=1{print $1"("a[$1]")",$2,$3,$4,$5,$6,$7}' Gene2Symbol.csv CountMatrix4DESeq1.csv > CountMatrix4DESeq1.symbol.csv 

cat CountMatrix4DESeq1.symbol.csv|awk 'BEGIN{FS=OFS=","}{print $1,$2,$3,$5,$6,$7}'  > temp

mv temp CountMatrix4DESeq1.symbol.csv 

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ../R_input/CountMatrix4DESeq1.symbol.csv ../R_input/colData2.csv EV ./


Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ../R_input/CountMatrix4DESeq1.symbol.csv ../R_input/colData.csv NC ./


#ensemble to symbol
Rscript GetGeneEntreID.R Homo_sapiens hsa


awk '
NR==FNR{a[$1]=$2;next}
{print $1,$2,a[$1] }' Gene2Symbol.txt gene_symbol.list > out

###转symbol

#第三列为合并。作为后期索引
cat Gene2Symbol.txt|awk '
BEGIN{FS=OFS="\t"}
{print $1,$2,$1"("$2")"}' >gene_symbol.list



#ensemble 对应 
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$3]=$2;next}
{print $1,$2,$3,a[$1]}' Homo_sapiens_GeneNames.txt gene_symbol.list |head

#第四列无则取第二列
cat out.list|awk '
BEGIN{FS=OFS="\t"}
$4==""&&$2!~/\./{print $1,$2,$3,$2}
$4!=""{print $0}
$4==""&&$2~/\./{print $0} ' >out.list2

#count 第四列有结果的
cat out.list2|awk '$4!=""{num++}END{print num}'
#40464

#第二列 不带点的
#39296
cat out.list|awk '
BEGIN{FS=OFS="\t"}
$2!~/\./{num++}
END{print num}'

##第四列为空的提出
cat out.list2|awk '
BEGIN{FS=OFS="\t"}
$4==""{print $0}'>no_symbol.list

##conversion 加入新列
for file in `ls */*`
do
	filename=${file/.txt/}
	awk '
	BEGIN{FS=OFS="\t"}
	NR==FNR{a[$3]=$4;next}
	FNR==1{print "Symbol\t"$0}
	FNR>1{print a[$1],$0}' ../../out.list2 $file > $filename.symbol.txt
done

##GO&KEGG 第9列拆分开

for file in `ls *.txt`
do
	filename=${file/_Symbol.txt/}

	cp $file temp1

	#GO 第五列格式有问题
	# cat $file|awk '
	# BEGIN{FS=OFS="\t"}
	# NR==1{print $0}
	# NR>1{split($5,x," ");
	# print $1,$2,$3,$4,x[1],x[2],x[3],$6,$7}' > temp1
 
	#KEGG 去掉第9列
	# cat $file|awk 'BEGIN{FS=OFS="\t"}
	# {print $1,$2,$3,$4,$5,$6,$7,$8,$10}' > temp1

	#拆分，放在最后
	cat temp1|awk '
	BEGIN{FS=OFS="\t"}
	{split($8,x,"/");
	for (i in x) print $0,x[i]}' > temp2

	#按第九列 对应symbol。放在最后
	awk '
	BEGIN{FS=OFS="\t"}
	NR==FNR{a[$3]=$4;next}
	FNR==1{print $0}
	FNR>1{print $0,a[$10]}' ../../../out.list2 temp2 > temp3

	#第二列相同则合并第11列
	cat temp3|awk 'BEGIN {FS=OFS="\t"}
	{
	    curr = $2
	    if (curr == prev) {
	        rec = rec "/" $11
	    }
	    else {
	        if (rec) print rec
	        rec = $0
	    }
	    prev = curr
	}
	END { if (rec) print rec }' > temp4

	cat temp4 |awk 'BEGIN{FS=OFS="\t"}
	NR==1{print "ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\tsymbolID"}
	NR>1{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$11}' > $filename.symbol.txt

done


file=KEGG_out_Symbol.txt


#heatmap 5 samples
awk '
BEGIN{FS=OFS=","}
NR==FNR{FS="\t";a[$3]=$4;next}
NR!=FNR && FNR==1 {FS=",";print $0}
NR!=FNR && FNR!=1 && a[$1]!=""{print a[$1],$2,$3,$4,$5,$6}' ../11.symbol/out.list2 CountMatrix4DESeq1.symbol.csv > CountMatrix4DESeq3.csv 
#duplicates
cat CountMatrix4DESeq3.csv|sort -t "," -k1,1 -u >temp
mv temp CountMatrix4DESeq3.csv
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ../R_input/CountMatrix4DESeq3.csv ../R_input/colData2.csv EV ./

#6 samples
awk '
BEGIN{FS=OFS=","}
NR==FNR{FS="\t";a[$3]=$4;next}
NR!=FNR && FNR==1 {FS=",";print $0}
NR!=FNR && FNR!=1 && a[$1]!=""{print a[$1],$2,$3,$4,$5,$6,$7}' out.list2 CountMatrix4DESeq1.symbol.csv > CountMatrix4DESeq3.csv 
#duplicates
cat CountMatrix4DESeq3.csv|sort -t "," -k1,1 -u >temp
mv temp CountMatrix4DESeq3.csv
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ../R_input/CountMatrix4DESeq3.csv ../R_input/colData.csv NC ./
