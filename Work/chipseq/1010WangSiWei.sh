'''
1.chipseq 
2.组蛋白修饰
3.新的参考基因组 建索引 （最好不从流程建立）
如果索引中中断，要把之前的临时文件都删除（fai，bt2，还有一个文件夹）。在重新建

'''


#rename WRXXXXX_R1.fastq.gz WRXXXXX_R2.fastq.gz

for smp in WR181364S WR181365S WR181366S WR181367S WR181368S WR181369S 
do 
	R1=`ls *${smp}*_R1*`
	R2=`ls *${smp}*_R2*`
	echo "sudo mv $R1 ${smp}_R1.fastq.gz"
	echo "sudo mv $R2 ${smp}_R2.fastq.gz"
done


#filter+QC

/data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_09/WangSiWei /data2/ClientData/2018_10/WangSiWei

#rename
while read -r old_name new_name; do
    rename $old_name $new_name *$old_name*.gz
done < ../filename


##conf prep
###gff2gtf


###chr_len.txt
seqkit fx2tab -l -g -n -i -H Mus_musculus.GRCm38.94.chr.fa |cut -f 1,4|sed '1d' > chr_len.txt

##kegg org
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/0Get_OrgDb_KEGG_Name.R "musculus"

"AH57974"
 kegg_code scientific_name common_name
14       mmu    Mus musculus       mouse

#-g GSIZE, --gsize GSIZE
Effective genome size. It can be 1.0e+9 or 1000000000,
or shortcuts:'hs' for human (2.7e9), 'mm' for mouse
(1.87e9), 'ce' for C. elegans (9e7) and 'dm' for
fruitfly (1.2e8), Default:hs

#coverge windown
5k  -  30M
50k - 300M 
10k    - 200M  

cat chr_len.txt|sort -k2,2nr|head  


#gtf 2 
1,4,5,9
cat Mus_musculus.GRCm38.94.chr.gtf|awk '
BEGIN {FS=OFS="\t"}
$3=="gene"{
	split($9,x,";")
	for (i in x) {
		if (x[i]~/gene_id/) {
			match(x[i],/\".*\"/,a)
			gsub(/"/,"",a[0])
			geneID=a[0];}
		if (x[i]~/gene_name/) {
			match(x[i],/\".*\"/,b)
			gsub(/"/,"",b[0])
			geneName=b[0];}
		}
	print $1,$4-3000,$5+3000,geneID"_"geneName

}' > MusGenesPlusMinus3kb.bed

##rpm
#wget http://hgdownload.soe.ucsc.edu/admin/exe/userApps.v350.src.tgz 
# conda install ucsc-wigtobigwig
# for smp in N N_input 
# bigWigToWig 

# bigWigToWig ../P_clean_sorted.bigWig P_clean_sorted.wig

#reads 和 peak 覆盖数统计 第10列（倒数第4列）
for smp in N P PH
do
	bedtools coverage -a ../${smp}_vs_${smp}_input_peaks.broadPeak -b ${smp}_clean.bed > ${smp}_peak_cov.txt
	bedtools coverage -a ../${smp}_vs_${smp}_input_peaks.broadPeak -b ${smp}_input_clean.bed > ${smp}_input_peak_cov.txt
done

#rpm=覆盖数/所有chr总reads数 * 1000000
#最后一列为rpm值
sed '1d' Statistic_Mapping.txt|cut -f 1,3 > temp
while read -r smp reads; do
    cat ${smp}_peak_cov.txt|cut -f 1,2,3,4,10|awk -v total="$reads" '
    BEGIN{FS=OFS="\t"}
    {rpm=$5/total*1000000;print $0,rpm}'|sed "1i chr\tstart\tend\tpeak_id\tnum\t${smp}_rpm" > ${smp}_peak_rpm.txt
done < temp
 
#overlap peaks between different group
mergePeaks -d 100 N_vs_N_input_peaks.broadPeak P_vs_P_input_peaks.broadPeak -prefix merge -venn P_vs_N_vennplot.pdf

mergePeaks -d 100 PH_vs_PH_input_peaks.broadPeak P_vs_P_input_peaks.broadPeak -prefix merge -venn PH_vs_P_vennplot.pdf

# #整合rpm值
# cat merge_N_vs_N_input_peaks.broadPeak_P_vs_P_input_peaks.broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > temp1
# #N
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$2]}' N_peak_rpm.txt temp1 > temp2
# #N_input
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$2]}' N_input_peak_rpm.txt temp2 > temp3

# #P
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$3]}' P_peak_rpm.txt temp3 > temp4
# #P_input
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$3]}' P_input_peak_rpm.txt temp4 > final.txt

# cat final.txt|awk 'BEGIN{FS=OFS="\t"}
# NR==1{print $0"\tN\tP"}
# NR>1{print $0,$4/$5,$6/$7}'|cut -f 1,8,9 > NvsP_rpm.txt



# cat merge_PH_vs_PH_input_peaks.broadPeak_P_vs_P_input_peaks.broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > temp1
# #PH
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$2]}' PH_peak_rpm.txt temp1 > temp2
# #PH_input
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$2]}' PH_input_peak_rpm.txt temp2 > temp3

# #P
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$3]}' P_peak_rpm.txt temp3 > temp4
# #P_input
# awk 'BEGIN{FS=OFS="\t"}
# NR==FNR{a[$4]=$6;next}
# {print $0,a[$3]}' P_input_peak_rpm.txt temp4 > final.txt

# cat final.txt|awk 'BEGIN{FS=OFS="\t"}
# NR==1{print $0"\tPH\tP"}
# NR>1{print $0,$4/$5,$6/$7}'|cut -f 1,8,9 > PHvsP_rpm.txt



#####fold enrichment 
cat merge_PH_vs_PH_input_peaks.broadPeak_P_vs_P_input_peaks.broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > PHvsP.txt
cat merge_N_vs_N_input_peaks.broadPeak_P_vs_P_input_peaks.broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > NvsP.txt

for smp in PH P N;do
	grep -v "#"  ${smp}_vs_${smp}_input_peaks.xls|cut -f 7,9|sed "1i $smp\tpeak_id" > ${smp}.enrich
done

#
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$2]}' PH.enrich PHvsP.txt > temp1

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$3]}' P.enrich temp1|cut -f 1,4,5 > PHvsP_Enrichment.txt

#
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$2]}' N.enrich NvsP.txt > temp1

awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$3]}' P.enrich temp1|cut -f 1,4,5 > NvsP_Enrichment.txt



####reads coverage
for smp in  N N_input P P_input PH PH_input 
do
	Rscript reads_coverage.R ${smp}_clean_plus_Count.txt ${smp}_clean_minus_Count.txt ${smp}_coverage.png
done



#峰 宽度分布
cd 2.CallPeak
for group in  N_vs_N_input P_vs_P_input PH_vs_PH_input 
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,4|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist.R temp 10000 img/${group}_peaks_length
done

#峰 富集倍数
for group in  N_vs_N_input P_vs_P_input PH_vs_PH_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,7|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 5 img/${group}_peaks_enrichment "Peak enrichment" "Count of peaks"
done

#峰 显著程度
for group in  N_vs_N_input P_vs_P_input PH_vs_PH_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,8|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 5 img/${group}_peaks_log "-log10 qvalue" "Count of peaks"
done









##比对统计 重新生成
for smp in N P PH
do
	echo "${smp},${smp}_clean_bowtie2.log" >> mapping_stat_list.txt
	echo "${smp}_input,${smp}_input_clean_bowtie2.log" >> mapping_stat_list.txt

	echo "${smp},${smp}_bam_stat.out" >> bam_stat_list.txt
	echo "${smp}_input,${smp}_input_bam_stat.out" >> bam_stat_list.txt
done

python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > Statistic_Mapping.txt

###Comm
bash /data1/script/deseq2+GO+KEGG/Rpipe/PlotEnrichHeatmapVenn.sh ./ PH P
bash /data1/script/deseq2+GO+KEGG/Rpipe/PlotEnrichHeatmapVenn.sh ./ P N





