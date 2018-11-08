#!/bin/bash
#cmscan结果整理和统计

#sRNA.sh cmscan.tblout /data1/Rfam/rfam_anno.txt ncRNA.stats ncRNA.length sRNA.gff "outData/sRNA/"

cat $2 | awk 'BEGIN {FS=OFS="\t"}{split($3,x,";");class=x[2];print $1,$2,$3,class}' > $6/rfam_anno_class.txt

grep -v "#" $1|awk '
BEGIN {FS=" "; OFS="\t"}
{if($20 != "=")print $2,$3,$4,$10,$11,$12,$17,$18}'|sed '1 i target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue' > $1.final.xls

#统计各RNA类别个数
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type=="") type="Others"; count[type]+=1;}END{for(type in count) print type, count[type];}' $6/rfam_anno_class.txt $1.final.xls > $3

#统计各RNA类别总长度
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type=="") type="Others"; if($6=="-")len[type]+=$4-$5+1;if($6=="+")len[type]+=$5-$4+1}END{for(type in len) print type, len[type];}' $6/rfam_anno_class.txt $1.final.xls > $4

#输出sRNA的全部记录
#加入第二列feature：sRNA
#调整至gff格式。并将start end顺序调整，end永远>start

awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type==" sRNA") print $0;}' $6/rfam_anno_class.txt $1.final.xls | awk 'BEGIN{FS=OFS="\t";}{print $3,"sRNA",$4,$5,$6,$1; }'|sed '1 i #Sequence\tFeature\tStart\tEnd\tStrand\tTarget_name'  |awk '
BEGIN{FS=OFS="\t"}
FNR!=1{
	if ($5=="+") {start=$3;end=$4}
	if ($5=="-") {start=$4;end=$3}
	print $1,"rfam",$2,start,end,"-",$5,$6}
FNR==1{print $1,"Source",$2,$3,$4,"Score",$5,$6}

' > $5 

