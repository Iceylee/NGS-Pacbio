#!/bin/bash

#0827 张毅无参项目 指定基因
#$1=gene_name $2=group_name

for gene in CP_G CP_GY GR_G GR_GY IR_G IR_GY OBP_G OR_G OR_GY SNM_G
do
	for group in FvsM LAvsFM SvsF
	do

		#根据第二列gene编号对应到nr的trinity编号，替换原来的第一列
		#去掉TRINITY_DN56740_c0_g1_i1的_i1
		awk 'BEGIN{FS=OFS="\t"}
		    NR==FNR{a[$2]=$1;next}
		    {print a[$2],$2,$3}' Coleoptera_nr.out Genes/${gene}.tsv |awk '
		    BEGIN{FS=OFS="\t"}
		    {split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4}' > temp1

		#提取各组的norm count数（按trinity编号）
		if [[ "$group" = "FvsM" || "$group" = "SvsF" ]]; then
		    awk '
		    BEGIN{FS=OFS="\t"}
		    NR==FNR{a[$1]=$2;next}
		    FNR==1{print "id\t"$0}
		    FNR>1{if ($1 in a) print a[$1],$0}' temp1  ${group}_norm-count-matrix.txt|cut -f 1,3,4,5,6 > temp2
		else
		    awk '
		    BEGIN{FS=OFS="\t"}
		    NR==FNR{a[$1]=$2;next}
		    FNR==1{print "id\t"$0}
		    FNR>1{if ($1 in a) print a[$1],$0}' temp1  ${group}_norm-count-matrix.txt|cut -f 1,3,4,5,6,7,8,9,10 > temp2
		fi

		#将Trinity编号换成gene name
		awk '
			BEGIN{FS=OFS="\t"}
			NR==FNR{a[$1]=$2;next}
			FNR==1{print $0}
			FNR>1{$1=a[$1];print $0}' Genes/gene_name.tab temp2 > temp3

		#name相同的加上序号 :2 来区分
		cat temp3|awk 'BEGIN{FS=OFS="\t"}
			{num[$1]=num[$1]+1;
			 if(num[$1]>=2)
					{idx = num[$1];
					$1=$1": "idx};
					print $0}' > PlotData/${group}_${gene}.tsv

		#作图
		Rscript heatmap.R PlotData/${group}_${gene}.tsv Heatmap/${group}_${gene}_ colData_${group}.csv

	done
done