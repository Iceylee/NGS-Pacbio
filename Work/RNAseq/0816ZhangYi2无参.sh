'''
1.追加4个样本。70 71 72 73
和之前的4个样本一起重新组装 
2.差异分析有三组比较
MS+FS vs FF+MF
FF+FS vs MF+MS
L1 L2 A1 A2 vs FF FS MF MS
3.增加分析：1）keggPathway 2）fasta注释 3）基因家族热图
4.KEGG，GO，Pathway，基因家族仅做鞘翅目（NR注释来过滤）
5.结果：差异表达结果在各自文件夹中 
'''

#QC

for smp in WR181070S WR181071S WR181072S WR181073S
do
	echo "mv *$smp*R1* ${smp}_R1.fastq.gz"
	echo "mv *$smp*R2* ${smp}_R2.fastq.gz"
done

bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_08/zhangyi  /data2/ClientData/2018_08/ZhangYi2

#conf
cat sample_list.txt|cut -f 1,2 > sample_exp_list.txt

cat sample_list.txt|cut -f 2|xargs -I {} echo {}/RSEM.genes.results > quant_files.txt


# gene_trans_map 获取 （流程跑完，发现该文件为空）
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/support_scripts/get_Trinity_gene_to_trans_map.pl /data2/ClientData/2018_08/ZhangYi2/1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa > /data2/ClientData/2018_08/ZhangYi2/1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa.gene_trans_map

bowtie2-build /data2/ClientData/2018_08/ZhangYi2/1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa /data2/ClientData/2018_08/ZhangYi2/1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa.bowtie2

#
cd /data2/ClientData/2018_08/ZhangYi2/2.Mapping/; for i in $(find . -name '*mapping.log');do echo `dirname $i|awk -F'/' '{print $2}'`,`basename $i` >> mapping_stat_list.txt ;done

cd /data2/ClientData/2018_08/ZhangYi2/2.Mapping/; python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > Statistic_Mapping.txt

#
cd %s; find . -name 'RSEM.genes.results' | sort | tee quant_files.txt


#按CountMatrix的样本顺序 调整 coldata的顺序




#分组差异分析
#1.SvsF
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ./CountMatrix4DESeq_SvsF.csv ./colData_SvsF.csv F /data2/ClientData/2018_08/ZhangYi2//4.DiffExprGene/DESeq2/SvsF/ TRUE

#2.FvsM
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ./CountMatrix4DESeq_FvsM.csv ./colData_FvsM.csv M /data2/ClientData/2018_08/ZhangYi2//4.DiffExprGene/DESeq2/FvsM/ TRUE

#3.LAvsFM
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ./CountMatrix4DESeq_LAvsFM.csv ./colData_LAvsFM.csv FM /data2/ClientData/2018_08/ZhangYi2//4.DiffExprGene/DESeq2/LAvsFM/ TRUE



##Gene Family

cat all.list|sed "s/\[.*\]//"|sed "s/PREDICTED: //"|sed "s/LOW QUALITY PROTEIN: //" > temp
cat temp|awk '
BEGIN{FS=" ";OFS="\t"}{$1=$1}1'| awk '
BEGIN{FS=OFS="\t"}
{printf $1"\t";for(i=2;i<=NF;i++)printf$i" ";print ""}' >gene_name.tab

cat gene_name.tab|sort -k1,1 -u > temp
mv temp gene_name.tab
#大写改成小写
#删除refName字样

#脚本gene-name.sh

#928 输出trinity编号而不是name
gene-name2.sh









##Statistics
for file in Annotation_kog.out Annotation_nr.out AllGene_GO_Annotation.txt AllGene_KEGG_Annotation.txt Annotation_uniprot.out  
do
	cat $file|sort -k1,1 -u|wc -l
done


for file in Annotation_kog.out Annotation_nr.out Annotation_uniprot.out GO_Annotation.txt
   56876 Annotation_kog.out
  102829 Annotation_nr.out
   62627 Annotation_uniprot.out
   46497 AllGene_KEGG_Annotation.txt
  54233 GO_Annotation.txt









##仅用鞘翅目基因3w
Annotation_nr.out #10w
cat Annotation_nr.out|grep "Aphis gossypii"|wc -l #80
cat Annotation_nr.out|grep "Acyrthosiphon pisum"|wc -l #12244
cat Annotation_nr.out|grep "Myzus persicae"|wc -l #1885
cat Annotation_nr.out|grep "Sitobion avenae"|wc -l #6
cat Annotation_nr.out|grep "Tribolium castaneum"|wc -l #12759

cat Annotation_nr.out|grep "Tribolium castaneum" > Tribolium_castaneum.nr.list

#提取第一列geneID和最后一列中括号中的物种
cat Annotation_nr.out|awk 'BEGIN{FS=OFS="\t"}{
	match($13,/\[.*\]/,a);
	gsub(/\[/,"",a[0]);
	gsub(/\]/,"",a[0]);
	print $1,a[0]}' > gene_species.tab

#取 Coleoptera目 36981 37%
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$1;next}
{if ($2 in a){print $1,$2}}' Coleoptera.txt gene_species.tab > Coleoptera_gene.tab

#Coleoptera nr 注释
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$1;next}
{if ($1 in a){print $0}}' Coleoptera_gene.tab Annotation_nr.out > Coleoptera_nr.out

#16803
Hemiptera

#2100
Lepidoptera

##鞘翅目重做KEGG和GO
#/data2/ClientData/2018_08/ZhangYi2/4.DiffExprGene/DESeq2/3.DiffExprGene 软链接到筛选鞘翅目的sig genes，然后重跑Enrichment脚本
group=LAvsFM
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{split($1,x,"_");
	$1=x[1]"_"x[2]"_"x[3]"_"x[4];
	a[$1]=$1;next}
{if ($1 in a)print $0}' Coleoptera_gene.tab ${group}_sig_genes_exprData.txt > ${group}_Coleoptera_sig_genes_exprData.txt

#鞘翅目重做基因家族热图
#将脚本gene-name.sh的15行的Annotation_nr.out改成Coleoptera_nr.out，重跑脚本。

#plot报错
Rscript /data1/script/GO_KEGG_Annotation/Plot_GO_KEGG.R SvsF_Coleoptera




#KEGG pathway

#取q最小的前20
group=FvsM
(head -n 1 ${group}_Coleoptera_KEGG_Enrichment.txt && tail -n +2 ${group}_Coleoptera_KEGG_Enrichment.txt|sort -t $'\t' -k5 -g) |head -26 > ../${group}_Coleoptera_KEGG_Enrichment.txt


#LAvsLM 19
#SvsF 21
#FvsM 8

#test
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/pathview.R ./4.DiffExprGene/DESeq2/3.DiffExprGene/ ./6.EnrichmentAnalysis/ ENSEMBL FALSE 


###3.fasta annotation
#nr.out 13
ln -s ../5.Annotation/Annotation_nr.out ./
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' Annotation_nr.out  Trinity.fasta >Trinity_anno.fasta

#real pep nr
pep=../7.TransDecoder/Trinity_CD-HIT_0.9.fa.transdecoder.pep
out=nr_pep.out
dbdir=/data1/DiamondDB/nr
log=nr_pep.log

echo "nohup /data1/software/diamond-linux64/diamond blastp -q $pep -d $dbdir -o ./$out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &"

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr_pep.out Trinity_CD-HIT_0.9.fa.transdecoder.pep >Trinity_anno.pep.fasta


#####
for smp in L1 L2 A1 A2 FF FS MF MS 
do
	ln -s ../../2.Mapping/${smp}/bowtie2.bam ${smp}_bowtie2.bam
done
