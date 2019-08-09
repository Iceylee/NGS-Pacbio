#!/bin/bash

if [ $# != 3 ];then
        echo
        echo "          USAGE: $0 CallPeakPath Sample1 Sample2"
        echo "          e.g.: $0 ./2.CallPeak A1 A2"
        echo
        exit 1;
fi

CallPeakPath=$1
smp1=$2
smp2=$3

if [  ! -d $CallPeakPath/CommDiffPeaks ]; then
	mkdir CommDiffPeaks
fi

cd $CallPeakPath/CommDiffPeaks

group=${smp1}vs${smp2}

ln -s ../${smp1}_vs_${smp1}_input_peaks.broadPeak ./
ln -s ../${smp2}_vs_${smp2}_input_peaks.broadPeak ./

#find overlap peaks between different group
mergePeaks -d 100 ${smp1}_vs_${smp1}_input_peaks.broadPeak ${smp2}_vs_${smp2}_input_peaks.broadPeak -prefix ${smp1}vs${smp2} -venn $group.venn.stats 

#heatmap
cat ${group}*broadPeak*broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > ${group}.txt

grep -v "#"  ../${smp1}_vs_${smp1}_input_peaks.xls|cut -f 7,9|sed "1i $smp1\tpeak_id" > ${smp1}.enrich
grep -v "#"  ../${smp2}_vs_${smp2}_input_peaks.xls|cut -f 7,9|sed "1i $smp2\tpeak_id" > ${smp2}.enrich

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$2]}' $smp1.enrich $group.txt > temp1

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$3]}' $smp2.enrich temp1|cut -f 1,4,5 > ${group}_Enrichment.txt

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/heatmap.R $group


#venn plot
smp1_stats=`cat $group.venn.stats|grep "$(printf '\t')${smp1}_vs_${smp1}_input_peaks.broadPeak$"|grep -E -o "$(printf '\t')[0-9]+$(printf '\t')"|grep -E -o "[0-9]+"`
smp2_stats=`cat $group.venn.stats|grep "$(printf '\t')${smp2}_vs_${smp2}_input_peaks.broadPeak$"|grep -E -o "$(printf '\t')[0-9]+$(printf '\t')"|grep -E -o "[0-9]+"`
comm_stats=`cat $group.venn.stats|grep "X$(printf '\t')X"|grep -E -o "$(printf '\t')[0-9]+$(printf '\t')"|grep -E -o "[0-9]+"`

echo "Rscript /data1/script/deseq2+GO+KEGG/Rpipe/vennplot.R $smp1_stats $smp2_stats $comm_stats $smp1 $smp2"

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/vennplot.R $smp1_stats $smp2_stats $comm_stats $smp1 $smp2


#delete useless files
rm -f ${smp1}.enrich ${smp2}.enrich temp1 ${group}.txt
