#!/bin/bash

#从cog注释结果统计cog分类，用于作图

# cog.sh /data1/COGDatabase/cog2003-2014.csv /data1/COGDatabase/cognames2003-2014.tab /data1/COGDatabase/fun2003-2014.tab cog.out  cog_count_fun.txt ./

#key第一列，需要第7列； 按照第二列的 gi|91777738|ref|YP_552946.1| 第二个查询
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$7;next}{split($2,x,"|");y=x[2];print $0, a[y]}' FS="," $1 FS="\t" $4 >$6/temp

#key第一列，需要第二列和第三列；按照14列COG编号查询
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{print $0, a[$14]}' $2 $6/temp > $6/temp2


#按照第15列分组计数。用于作图 (第二列，按数字，从大到小倒序排列)
cat $6/temp2 | awk 'BEGIN{FS=OFS="\t"}{count[$15]+=1}END{for(id in count)print id,count[id]}'|sort -k2nr,2 >$6/count.id
#多个类的拆分开，并加入各类的计数
cat $6/count.id | awk '
	BEGIN{OFS="\t"}
	{
		split($1,x,"")
		for (i in x){
			class = x[i]
			count[class]+=$2
		}
	}
	END{for (id in count)print id,count[id]}' |sort -k1,1 > $6/count.split.id

#count id加上一列，说明分类的意义；key第一列，需要第二列；按照第一列查询
awk 'BEGIN{OFS="\t"}FNR==NR{a[$1]=$2;next}{print $0, a[$1]}' FS="\t" $3 FS="\t" $6/count.split.id > $5

rm -f $6/temp $6/temp2 $6/count.split.id $6/count.id
