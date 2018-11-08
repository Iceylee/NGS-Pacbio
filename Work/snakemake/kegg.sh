#!/bin/bash

#从kaas结果提取kegg注释，并统计用于作图
# kegg.sh q00000.keg kegg.anno.out kegg_count.out

keg=$1
anno=$2
count=$3


cat $keg |grep "^D"|sed 's/^D      //g'|awk 'BEGIN{FS="; ";OFS="\t"}{sub(" ",":",$2);print $1,$2";"$3}' |sort -k1,1 |awk '
BEGIN{FS=OFS="\t"}{split($2,x,":");print $1,$2,x[1]}'|awk 'BEGIN{FS=OFS="\t"}!a[$1,$3]++'|cut -f 1,2 > $anno



cat $keg |awk 'BEGIN{FS=OFS="\t"}
/^A/{A_name=$0}
/^B  /{B_name=$0}
/^D/{D_name=$0;print A_name,B_name,D_name}' |awk '{FS=OFS="\t"}{split($3,x,";");split(x[1],y,"_");print $1,$2,$3,y[3]}' |sort -k3 -u |awk 'BEGIN{FS=OFS="\t"}{a[$2]++;b[$2]=$1}END{for (i in a)print b[i],i,a[i]}' |sed -r 's/A[0-9]+ //g'|sed -r 's/B  [0-9]+ //g'| sort -k3,3 -n > $count