#!/bin/bash

#按照tab，将第二列转换回第一列

$1 $2

if [  -f $1/$2 ]; then

	filename=${file/_ID_type.tab/}
	tabFile=${filename}_ID_type.tab
	keggFile=${filename}_KEGG_out.txt
	awk '
	BEGIN{FS=OFS="\t"}
	NR==FNR{a[$2]=$1;next}
	{	str=""
		split($8,x,"/");
		for (i in x) {
			p = x[i]
			if(str!=""){str = str"/"a[p]}
			if(str==""){str = a[p]}
			}
		print $1,$2,$3,$4,$5,$6,$7,str,$9}' $tabFile $keggFile | sed "1 d"|sed "1 i KEGG_ID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount" > temp

	mv -f temp $keggFile
	rm -f $tabFile
	
fi