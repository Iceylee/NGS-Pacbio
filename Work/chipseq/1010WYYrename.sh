

for smp in A1 A2 A3
do
	ln -s ../bk/1.Mapping/${smp}_clean_sorted.bam ${smp}_input_clean_sorted.bam
	ln -s ../bk/1.Mapping/${smp}_input_clean_sorted.bam ${smp}_clean_sorted.bam
done

for smp in A1 A2 A3
do
	ln -s ../../bk/1.Mapping/Depth_Count_50k/png/${smp}_coverage.png ${smp}_input_coverage.png
	ln -s ../../bk/1.Mapping/Depth_Count_50k/png/${smp}_input_coverage.png ${smp}_coverage.png
done

for smp in A1 A2 A3
do
	ln -s ../bk/1.Mapping/${smp}_clean_sorted.bam.bai ${smp}_input_clean_sorted.bam.bai
	ln -s ../bk/1.Mapping/${smp}_input_clean_sorted.bam.bai ${smp}_clean_sorted.bam.bai

done


#diffReps
#rename
for smp in A1 A2 A3
do
	mv ${smp}_input.bed ${smp}.bed
done

nohup perl /data1/software/diffreps/bin/diffReps.pl -tr A1.bed -co A2.bed -re A1vsA2_diff.txt -me gt -ch chr_len.txt &

nohup perl /data1/software/diffreps/bin/diffReps.pl -tr A1.bed -co A3.bed -re A1vsA3_diff.txt -me gt -ch chr_len.txt &

nohup perl /data1/software/diffreps/bin/diffReps.pl -tr A2.bed -co A3.bed -re A2vsA3_diff.txt -me gt -ch chr_len.txt &

#注释差异peak
HumanGenesPlusMinus3kb=/data1/GenomicDatabases/Human/Ensembl/HumanGenesPlusMinus3kb.bed
for group in A1vsA3 A1vsA2 A2vsA3
do
	intersectBed -wa -a $HumanGenesPlusMinus3kb -b ${group}_diff.txt.hotspot > diff_${group}_Genesat3KborlessfromPeaks.txt

	cat diff_${group}_Genesat3KborlessfromPeaks.txt|cut -f 4 |awk '{split($1,x,"_");print x[1]}'|sort -k1,1 -u|sed '1 i GeneID'| > ${group}_sig_genes_exprData.txt
done

#KEGG GO
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa ENSEMBL ENTREZID /data2/ClientData/2018_09/WengYuanYuan/6.DiffPeak/diffReps/

#diffMotifs
GenomeFa=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.fa
for group in A1vsA3 A1vsA2 A2vsA3
do
	findMotifsGenome.pl ${group}_diff.txt.hotspot $GenomeFa ${group} -len 8,10,12
done


##frag size
Rscript A1_vs_A1_input_model.r
Rscript A2_vs_A2_input_model.r
Rscript A3_vs_A3_input_model.r


##3个分布图
##peak len
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,4|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist.R temp 10000 img/${group}_peaks_length
done

#峰 富集倍数
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,7|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 10 img/${group}_peaks_enrichment "Peak enrichment" "Count of peaks"
done

#峰 显著程度
for group in  A1_vs_A1_input A2_vs_A2_input A3_vs_A3_input
do
	cat ${group}_peaks.xls|grep -v "#"|cut -f 9,8|awk 'BEGIN{FS=OFS="\t"}{print $2,$1}'|sed '1,2d'  > temp

	Rscript /data1/script/deseq2+GO+KEGG/Rpipe/length_hist_5_argu.R temp 10 img/${group}_peaks_log "-log10 qvalue" "Count of peaks"
done


#snp rename
#先将input都存到temp中
#rename
while read -r old_name new_name; do
    rename $old_name $new_name $old_name*
done < filename

#IGV
mapping中的bigwig文件，IP和input同时导入

