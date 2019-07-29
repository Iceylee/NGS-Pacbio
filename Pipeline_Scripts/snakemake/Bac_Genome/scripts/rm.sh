#!/bin/bash

sed '1,3 d' $1|awk 'BEGIN{FS=" ";OFS="\t"}{locus=$5"-"$6"-"$7;print locus,$11}' > temp1
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}!/#/{locus=$1"-"$4"-"$5;print $0,a[locus]}' temp1 $1.gff |grep -v "#" |awk 'NR==1{print "Sequence\tSource\tMethod\tStart\tEnd\tScore\tStrand\tPhase\tTarget\tRepeat_class"}NR>1{print $0}' > $2