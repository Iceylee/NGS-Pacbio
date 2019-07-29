'''
1. 6样本，有一个样本用的补测的数据。
2. 样本名称不要有下划线。会导致比对结果统计文件出问题。
3. 组蛋白修饰 H3K18cr抗体
4. wrongName:IP与input搞颠倒了。bam开始重跑流程。cleanData文件（按原来的错误映射）重新更改到公司编号。SNP和coverage图换成正确的名字。（注意软链接的指向）

'''


#rename WRXXXXX_R1.fastq.gz WRXXXXX_R2.fastq.gz

for smp in WR181358S WR181359S WR181360S WR181361S WR181362S WR181363S
do 
	R1=`ls *${smp}*_R1*`
	R2=`ls *${smp}*_R2*`
	echo "mv $R1 ${smp}_R1.fastq.gz"
	echo "mv $R2 ${smp}_R2.fastq.gz"
done

mv S035_whbioacme-A_WR181363S_BH5VCVDSXX_S7_L002_R1_001.fastq.gz WR181363S_new_R1.fastq.gz
mv S035_whbioacme-A_WR181363S_BH5VCVDSXX_S7_L002_R2_001.fastq.gz WR181363S_new_R2.fastq.gz

#filter+QC

/data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_09/Wenyuanyuan /data2/ClientData/2018_09/WengYuanYuan

Step1_FastQC_analysis_PE.sh /Databackup/2018_09/Wenyuanyuan /data2/ClientData/2018_09/WengYuanYuan

for i in WR181358S WR181359S WR181360S WR181361S WR181362S WR181363S
do
	/data1/software/Reseqtools/iTools_Code/iTools Fqtools fqcheck -InFq1 clean_data_temp/${i}_clean_R1.fastq.gz -OutStat1 FastQC/clean/${i}_clean_R1 -InFq2 clean_data_temp/${i}_clean_R2.fastq.gz -OutStat2 FastQC/clean/${i}_clean_R2
done

for i in WR181363S_new
do
	/data1/software/Reseqtools/iTools_Code/iTools Fqtools fqcheck -InFq1 clean_data_temp/${i}_clean_R1.fastq.gz -OutStat1 FastQC/clean/${i}_clean_R1 -InFq2 clean_data_temp/${i}_clean_R2.fastq.gz -OutStat2 FastQC/clean/${i}_clean_R2
done

#63
/data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_09/Wenyuanyuan/63 /Databackup/2018_09/Wenyuanyuan/63_clean_fastqc

#rename
while read -r old_name new_name; do
    rename $old_name $new_name *$old_name*.gz
done < ../filename

#rename _h _half *.png
#rename 's/$old_name/$new_name/' *$old_name*.txt

# pipeline
#流程需要用的文件
/data1/script/Pipeline/ChIP-Seq  #conf,py,PeakAnnotation.R,SNP_Indels.py

#之前测试过的目录
/data2/ClientData/2018_06/HuangQiLai/Chip-Seq_PipelineTest


#下划线命名样本 导致 比对结果为空  
#重跑
for smp in A1 A1_input A2 A2_input A3 A3_input
do

	echo "/data1/software/RSeQC-2.6.4/scripts/bam_stat.py -i /data2/ClientData/2018_09/WengYuanYuan/1.Mapping/${smp}_clean_sorted.bam > /data2/ClientData/2018_09/WengYuanYuan/1.Mapping/${smp}_bam_stat.out"
done > cmd

bash cmd


python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > Statistic_Mapping.txt








##Motif


Contrast = j + "_vs_" + i
NarrowPeak = OutPut + "2.CallPeak/" + Contrast + "_peaks.narrowPeak"
NarrowPeak_fa = OutPut + "5.FindMotif/" + Contrast + "_peaks.fa"

Command11 = "bedtools getfasta -fi %s -bed %s -name > %s" % (GenomeFa, NarrowPeak, NarrowPeak_fa)
print "--** Get Peak Fa **--Command11: \n%s\n\n" % Command11
os.system(Command11)

NarrowPeak_homer_bed = OutPut + "5.FindMotif/" + Contrast + "_peaks_homer.bed"
MotifDir = OutPut + "5.FindMotif/"
Command12_1 = "NarrowPeak2HomerBed.sh %s > %s" % (NarrowPeak, NarrowPeak_homer_bed)
#awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $1


Command12_2 = "findMotifsGenome.pl %s %s %s -len 8,10,12" % (NarrowPeak_homer_bed, GenomeFa, MotifDir)
print "--** Get Motif **--Command12_1: \n%s\n\n" % Command12_1
print "--** Get Motif **--Command12_2: \n%s\n\n" % Command12_2

##call peak
###predict fragment size
macs2 predictd -i A2.bed -g hs -m 10 30

for smp in A1 A1_input A2 A2_input A3 A3_input 
do
	echo "macs2 predictd -i ${smp}.bed -g hs -m 5 50 &> $smp.pre &"
done

for smp in A1 A1_input A2 A2_input A3 A3_input
do
	cat $smp.pre|grep "predicted fragment length"
done



Command9 = "macs2 callpeak -t %s -c %s -m 5 50 -p 1e-5 -f BAM -g %s -n %s --outdir %s -B --broad --nomodel --extsize 236 2>%s"

#frament size plot
Rscript A3_vs_A3_input_model.r

##coverage plot

#split strand and output depth
for smp in A1 A1_input A2 A2_input A3 A3_input 
do
	# foward
	echo "samtools view -F 16 -h ${smp}_clean_sorted.bam |samtools sort -@ 15 - > ${smp}_plus.bam"

	# reverse 
	echo "samtools view -f 16 -h ${smp}_clean_sorted.bam |samtools sort -@ 15 - > ${smp}_minus.bam"

	echo "samtools depth ${smp}_minus.bam > ${smp}_minus_depth.txt"

	echo "samtools depth ${smp}_plus.bam > ${smp}_plus_depth.txt"

done > depth.cmd

bash depth.cmd &> depth.log

# genomeCoverageBed -ibam sample.bam -d -g genome.bed

# bedtools genomecov -ibam LM_sorted.bam -d >LM.dep_base

# samtools mpileup

#R plot 输入数据：1k窗口median
#Log2


##macs diffbdg
egrep "tags after filtering in treatment|tags after filtering in control" A1_vs_A1_input_peaks.xls
# tags after filtering in treatment: 31589830
# tags after filtering in control: 43607487
egrep "tags after filtering in treatment|tags after filtering in control" A2_vs_A2_input_peaks.xls
# tags after filtering in treatment: 28160519
# tags after filtering in control: 30951484
egrep "tags after filtering in treatment|tags after filtering in control" A3_vs_A3_input_peaks.xls
# tags after filtering in treatment: 29105983
# tags after filtering in control: 27813725

#只看control
# -g 60 -l 120  #好像要和callpeak的extsize一致？？ 默认是100和200
macs2 bdgdiff --t1 A1_vs_A1_input_treat_pileup.bdg --c1 A1_vs_A1_input_control_lambda.bdg --t2 A2_vs_A2_input_treat_pileup.bdg --c2 A2_vs_A2_input_control_lambda.bdg --d1 43607487 --d2 30951484 --o-prefix diff_A1_vs_A2

macs2 bdgdiff --t1 A1_vs_A1_input_treat_pileup.bdg --c1 A1_vs_A1_input_control_lambda.bdg --t2 A3_vs_A3_input_treat_pileup.bdg --c2 A3_vs_A3_input_control_lambda.bdg --d1 43607487 --d2 27813725 --o-prefix diff_A1_vs_A3

macs2 bdgdiff --t1 A2_vs_A2_input_treat_pileup.bdg --c1 A2_vs_A2_input_control_lambda.bdg --t2 A3_vs_A3_input_treat_pileup.bdg --c2 A3_vs_A3_input_control_lambda.bdg --d1 30951484 --d2 27813725 --o-prefix diff_A2_vs_A3

#diff gene anno
#合并cond1和cond2
for group in A1_vs_A3 A1_vs_A2 A2_vs_A3
do 
	sed '1 d' diff_${group}_c3.0_cond1.bed |awk 'BEGIN{FS=OFS="\t"}{
		$6="cond1";print $0}' > diff_${group}.bed
	sed '1 d' diff_${group}_c3.0_cond2.bed |awk 'BEGIN{FS=OFS="\t"}{
		$6="cond2";print $0}' >> diff_${group}.bed
done

##方法一 ：这个注释的是转录本 ENST。。。不能用于后续的GO和KEGG注释
Genomebed=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.bed
GenomeGTF=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf
for group in A1_vs_A3 A1_vs_A2 A2_vs_A3
do
	annotatePeaks.pl diff_${group}.bed $Genomebed -gtf $GenomeGTF > ${group}_peaks_annotation.txt
done

for group in A1_vs_A3 A1_vs_A2 A2_vs_A3
do
	sed '1d' ${group}_peaks_annotation.txt|cut -f 11 |sort -k1,1 -u|sed '1 i GeneID'|sed '$ d'> ${group}_sig_genes_exprData.txt
done

##方法二：
HumanGenesPlusMinus3kb=/data1/GenomicDatabases/Human/Ensembl/HumanGenesPlusMinus3kb.bed
for group in A1_vs_A3 A1_vs_A2 A2_vs_A3
do
	intersectBed -wa -a $HumanGenesPlusMinus3kb -b diff_${group}.bed > diff_${group}_Genesat3KborlessfromPeaks.txt
done

for group in A1_vs_A3 A1_vs_A2 A2_vs_A3
do
	cat diff_${group}_Genesat3KborlessfromPeaks.txt|cut -f 4 |awk '{split($1,x,"_");print x[1]}'|sort -k1,1 -u|sed '1 i GeneID'| > ${group}_sig_genes_exprData.txt
done


##GO&KEGG

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa ENSEMBL ENTREZID /data2/ClientData/2018_09/WengYuanYuan/6.DiffPeak/GOKEGG/





##diffReps 针对无重复的peak
#treatment:tr实验 control ：co对照
# (nb=Negative binomial; gt=G-test; tt=T-test; cs=Chi-square test).

#bam2bed
for smp in A1 A1_input A2 A2_input A3 A3_input 
do
	echo "bedtools bamtobed -i ../1.Mapping_bk/${smp}_clean_sorted.bam > ${smp}.bed"
done >cmd

nohup cmd &> log

#diffreps
perl /data1/software/diffreps/bin/diffReps.pl -tr A1.bed -co A2.bed -re A1vsA2_diff.txt -me gt -ch chr_len.txt

nohup perl /data1/software/diffreps/bin/diffReps.pl -tr A1.bed -co A3.bed -re A1vsA3_diff.txt -me gt -ch chr_len.txt &> 13.log &

perl /data1/software/diffreps/bin/diffReps.pl -tr A2.bed -co A3.bed -re A2vsA3_diff.txt -me gt -ch chr_len.txt &> 23.log &

#注释差异peak
HumanGenesPlusMinus3kb=/data1/GenomicDatabases/Human/Ensembl/HumanGenesPlusMinus3kb.bed
for group in A1vsA3 A1vsA2 A2vsA3
do
	intersectBed -wa -a $HumanGenesPlusMinus3kb -b ${group}_diff.txt.hotspot > diff_${group}_Genesat3KborlessfromPeaks.txt
done

#
for group in A1vsA3 A1vsA2 A2vsA3
do
	cat diff_${group}_Genesat3KborlessfromPeaks.txt|cut -f 4 |awk '{split($1,x,"_");print x[1]}'|sort -k1,1 -u|sed '1 i GeneID'| > ${group}_sig_genes_exprData.txt
done

#KEGG GO
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa ENSEMBL ENTREZID /data2/ClientData/2018_09/WengYuanYuan/6.DiffPeak/diffReps/

##diffMotif
GenomeFa=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.fa
findMotifsGenome.pl A1vsA2_diff.txt.hotspot $GenomeFa ./ -len 8,10,12




#峰 宽度分布
cd 2.CallPeak
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,4|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > ${group}_peaks_length.txt

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist.R ${group}_peaks_length.txt 10000 ${group}_peaks_length
done

#峰 富集倍数
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,7|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 5 img/${group}_peaks_enrichment "Peak enrichment" "Count of peaks"
done

#峰 显著程度
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,8|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 5 img/${group}_peaks_log "-log10 qvalue" "Count of peaks"
done


5k  -  30M
50k - 300M ??



#coverage plot
Rscript reads_coverage.R A3_plus_Count.txt A3_minus_Count.txt A3_coverage.png

for smp in  A1 A1_input A2 A2_input A3 A3_input 
do
	Rscript reads_coverage.R ${smp}_plus_Count.txt ${smp}_minus_Count.txt ${smp}_coverage.png
done

#UCSC
##染色体名称改为chr1 并只包含chr1-22 X Y
cat ../A1_vs_A1_input_peaks.bed | awk 'BEGIN{FS=OFS="\t"}
{
	$1="chr"$1
	for (i=1;i<=22;i++) {
		if($1=="chr"i){print $0}
	}
	if($1=="chrX"){print $0}
	if($1=="chrY"){print $0}
}' > A1_vs_A1_input_peaks4UCSC.bed


##rpm heatmap
cd 2.CallPeak/rpm2
#A1vsA3 A1vsA2 A2vsA3

#overlap peaks between different group
mergePeaks -d 100 A1_vs_A1_input_peaks.broadPeak A2_vs_A2_input_peaks.broadPeak -prefix A1vsA2

mergePeaks -d 100 A1_vs_A1_input_peaks.broadPeak A3_vs_A3_input_peaks.broadPeak -prefix A1vsA3

mergePeaks -d 100 A2_vs_A2_input_peaks.broadPeak A3_vs_A3_input_peaks.broadPeak -prefix A2vsA3

#####fold enrichment 
for group in A1vsA3 A1vsA2 A2vsA3;do
	cat ${group}*broadPeak*broadPeak|cut -f 1,9,10|sed '1d'|sed '1i merge_id\tpeak_id\tpeak_id' > ${group}.txt
done

for smp in A1 A2 A3;do
	grep -v "#"  ${smp}_vs_${smp}_input_peaks.xls|cut -f 7,9|sed "1i $smp\tpeak_id" > ${smp}.enrich
done

###
while read -r smp1 smp2 group; do
    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{a[$2]=$1;next}
    {print $0,a[$2]}' $smp1.enrich $group.txt > temp1

    awk 'BEGIN{FS=OFS="\t"}
    NR==FNR{a[$2]=$1;next}
    {print $0,a[$3]}' $smp2.enrich temp1|cut -f 1,4,5 > ${group}_Enrichment.txt
done < filename


##
for group in A1vsA3 A1vsA2 A2vsA3;do
	Rscript heatmap.R $group
done









