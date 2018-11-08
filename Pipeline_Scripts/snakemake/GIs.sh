#!/bin/bash

#GIs.sh 0 assembly.fasta gff


for file in `ls *_GIs.txt`
do
	name=${file/_GIs.txt/}
cat ${file} | awk -v contig="$name" 'BEGIN{OFS="\t"}{print "unitig_"contig,$0}' >> all.GIs.out
done


cat all.GIs.out|awk 'BEGIN{OFS="\t"}{$2="Genomic Island";print $0}' |sed "1 i Sequence\tFeature\tStart\tEnd" > islands.txt

# gff提取gene坐标
ln -s ../prodigal/Thermoactinomyces_sp.CDF.gff gff
awk 'BEGIN{OFS="\t"}
NR==FNR&&FNR>1{i++;sequence[i]=$1;start[i]=$3;end[i]=$4;next}
{for (i in start) 
	{if($4>=start[i]&&$5<=end[i]) print sequence[i],"Dimob","Gene",$4,$5,".",$7,$9,start[i],end[i],end[i]-start[i]+1}}' islands.txt gff> temp
cat temp | awk 'BEGIN {OFS="\t"}
	{print $0}'|sed "1 i #Sequence\tSource\tFeature\tGene_Start\tGene_End\tScore\tStrand\tGene_ID\tIsland_Start\tIsland_End\tIsland_Length" > GeneIslands.txt