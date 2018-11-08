'''
1.有参转录组，3个样本。
2.无生物学重复，用edgeR做差异分析
3.非模式生物，用eggnog做kegg和go注释
4.参考基因组的contig数高达130w（每条序列是一个基因），GTAK无法做SNP。用samtools做过滤和提取。


'''



#sample_list.txt


#exon_len.txt
python /data1/script/ExtraExonLenFromGTF.py GCA_000516895.1_LocustGenomeV1_genomic.gtf > Gene_ExonLen.txt


#Gene2Symbol.txt



#tab
LncRNAID_File = /data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.lncRNA.ID

#
mkdir R_input
colData.csv

screen -R LiDaQi
screen -ls
screen -r LiDaQi #重新进入
#ctr+A+D 退出
python TransSeq_WithGenome.py TransSeq_WithGenome.conf

#htseq 重跑 
## gtf 必须是gtf文件
gtf=/data3/ClientData/2018_07/LiDaQi/Genome/GCA_000516895.1_LocustGenomeV1_genomic.gtf
nohup htseq-count -f bam WR180755S_clean_sorted.bam $gtf -q > WR180755S.htseq.out &
nohup htseq-count -f bam WR180756S_clean_sorted.bam $gtf -q > WR180756S.htseq.out &
nohup htseq-count -f bam WR180757S_clean_sorted.bam $gtf -q > WR180757S.htseq.out &


Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_input/CountMatrix4DESeq.csv R_input/colData.csv Ck-24 ./ FALSE






#snpeff 构建基因组
#参数 Locusta_migratoria 为构建的基因组

cd /data1/software/snpEff/data
mkdir Locusta_migratoria #genes.gtf  sequences.fa  snpEffectPredictor.bin
cp /data3/ClientData/2018_07/LiDaQi/Genome/*.fa sequences.fa
cp /data3/ClientData/2018_07/LiDaQi/Genome/*.gtf genes.gtf

vim snpEff.config
#add following in """
 """
# Maize, version GCF_000005005.2_B73_RefGen_v4_genomic
Maize.genome : Maize
"""

#在snpEff目录
java -jar snpEff.jar build -gff3 -v  Locusta_migratoria2 #gff
java -jar snpEff.jar build -gtf22 -v Locusta_migratoria #gtf


python SNP_Indels.py Genome/GCA_000516895.1_LocustGenomeV1_genomic.fa ./ Locusta_migratoria







#gfold 差异分析
#两个实验组 1个对照：55
#gfold count -ann hg19_Ref.gtf -tag sample1.sam -o sample1.read_cnt


#筛选，GFOLD values等于0，表示该基因不差异；大于0，表示该基因上调；小于0，表示该基因下调。


samtools view WR180755S_clean_sorted.bam | gfold count -ann ../../Genome/GCA_000516895.1_LocustGenomeV1_genomic.gtf -tag stdin -o WR180755S.read_cnt
samtools view WR180756S_clean_sorted.bam | gfold count -ann ../../Genome/GCA_000516895.1_LocustGenomeV1_genomic.gtf -tag stdin -o WR180756S.read_cnt
samtools view WR180757S_clean_sorted.bam | gfold count -ann ../../Genome/GCA_000516895.1_LocustGenomeV1_genomic.gtf -tag stdin -o WR180757S.read_cnt

#sample2/sample1 s1为对照

gfold diff -s1 WR180755S -s2 WR180756S -suf .read_cnt -sc 0.01 -o WR180756S-VS-WR180755S.diff
gfold diff -s1 WR180755S -s2 WR180757S -suf .read_cnt -sc 0.01 -o WR180757S-VS-WR180755S.diff
#gfold 阈值 绝对值大于1
cat WR180756S-VS-WR180755S.diff|grep -v "#"|awk '$3>=1||$3<=-1{print $0}'>WR180756S-VS-WR180755S.sig.diff
cat WR180757S-VS-WR180755S.diff|grep -v "#"|awk '$3>=1||$3<=-1{print $0}'>WR180756S-VS-WR180755S.sig.diff


#RPKM
/data1/script/AnalysisRPKM.py Genome/Gene_ExonLen.txt sample_list.txt > 2.GenesExpress/AllSamplesRPKMValue.txt











##SNP
#### 1.Variant filtering ####
smp=Ck-24

variant_name=${smp}.samtools.raw.vcf
filter_name=${smp}.flt1.vcf
filter_log_name=${smp}.flt1.log
GenomeFa=/data3/ClientData/2018_07/LiDaQi/Genome/GCA_000516895.1_LocustGenomeV1_genomic.fa
GenomeDict=/data3/ClientData/2018_07/LiDaQi/Genome/GCA_000516895.1_LocustGenomeV1_genomic.dict

#java -jar /data1/software/picard/picard.jar CreateSequenceDictionary R=$GenomeFa O=$GenomeDict #1min

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T VariantFiltration  -R $GenomeFa -V $variant_name -window 35 -cluster 3 -filterName FSS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o $filter_name 2>$filter_log_name #3h

####1.1 samtools vcf filtering
#1.SNP的reads支持数不低于4；2.SNP的质量值(MQ)不低于20
#第8列的 DP4:4个值：1 比对结果和正链一致的reads数、2 比对结果和负链一致的reads数、3 比对结果在正链的variant上的reads数、4 比对结果在负链的variant上的reads数。可以设定 （value3 + value4）大于某一阈值，才算是variant
#第8列的MQ 

 
#DP 和 MQ 过滤
vcfutils.pl varFilter -Q 20 -d 4 $vcf > 

#提取
for smp in Ck-24 G-24 Ten-24
do
	#bcftools view -v snps ${smp}.flt1.vcf > ${smp}.SNPs.vcf
	#bcftools view -v indels ${smp}.flt1.vcf > ${smp}.INDELs.vcf

	java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s ${smp}.SnpEff.html -c /data1/software/snpEff/snpEff.config -v -ud 500 Locusta_migratoria ${smp}.flt1.vcf > ${smp}.SnpEff.vcf
done



			










 ####2. SNP/INDEL ####
SNP_name=${smp}.SNPs.vcf
INDEL_name=${smp}.INDELs.vcf

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType SNP -o $SNP_name # 1h

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType INDEL -o $INDEL_name #1h

 ###3.Annotation ###
Annota_Html=${smp}.SnpEff.html
Annota_vcf=${smp}.SnpEff.vcf
snpEff_SpeciesName=Locusta_migratoria

java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s $Annota_Html -c /data1/software/snpEff/snpEff.config -v -ud 500 $snpEff_SpeciesName $filter_name > $Annota_vcf












##################
##GO annotation & enrichment
##################

#eggnog注释: (GO KEGG COG/KOG)  (2.8W条序列--1hour)
emapper.py -i $pep --output eggnog.out -m diamond --cpu 16 # 6.27 9am- 6.29 1am 5w

emapper.py -i $nuc_fa --translate --output eggnog.out -m diamond --cpu 16

python /data1/software/eggnog-mapper/emapper.py -i GCA_000516895.1_LocustGenomeV1_genomic.fa --translate --output eggnog.out -m diamond --cpu 16

# 功能关联: (可跳过)
# emapper-ann -d GO eggnog.out.emapper.annotations > Name_GO_list.txt
# emapper-ann -d KEGG eggnog.out.emapper.annotations > Name_KEGG_list.txt
# emapper-eggnog -c COG eggnog.out.emapper.annotations > Name_COG_list.txt

#backgroud

python /data1/script/GO_KEGG_Annotation/GetGOID_Annotation.py eggnog.out.emapper.annotations /data1/NCBI/gene2go > GO_annotation.txt
   #tr名称去掉最后一个下划线后面内容，转为gene名称
cat GO_annotation.txt|awk '
BEGIN{FS=OFS="\t"}
{split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4}' > GO_annotation_gene.txt

#sig genes
cat ../../5.Diff/3.DiffExprGene/LvsA_sig_genes_exprData.txt|cut -f 1 > sig.list

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$1;next}
{if($1 in a) print $0}' sig.list GO_annotation_gene.txt > LvsA_GO_annotation.txt

#enrichment

python /data1/script/GO_KEGG_Annotation/Pyher_GOEnrichment_Analysis.py LvsA_GO_annotation.txt GO_annotation_gene.txt LvsA




###差异基因的KEGG注释
#提取差异基因的fasta序列
#id0 AVCP010000001.1
for group in 10-24vsCk-24 G-24vsCk-24
do
	fasta=GCA_000516895.1_LocustGenomeV1_genomic.fa
	sed '1d' ${group}_sig_genes_exprData.txt |awk '{match($1,/[0-9]+/,a);num=10000001+a[0];print "AVCP0"num".1"}' >${group}_sig.list
	cat $fasta|seqkit grep -f ${group}_sig.list > ${group}_sig.fasta
done


###差异基因的GO注释
group=G-24vsCk-24 #10-24vsCk-24
fasta=../${group}_sig.fasta

out=${group}.uniprot.out
dbdir=/data1/Uniprot/uniprot_sprot
log=uniprot.log

echo "nohup /data1/software/diamond-linux64/diamond blastx -q $fasta -d $dbdir -o ./$out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &"

#--sensitive 



##eggnog
#group=G-24vsCk-24 
#group=10-24vsCk-24
group=G-24vs10-24

echo "python /data1/software/eggnog-mapper/emapper.py -i ../${group}_sig.fasta --translate --output ${group}.eggnog.out -m diamond --cpu 16"

python /data1/script/GO_KEGG_Annotation/GetKEGGID_Annotation.py ${group}.eggnog.out.emapper.annotations /data1/KEGG/ko00000_kegg.2018.5.30.txt > ${group}.KEGG_Annotation.txt


cat *.keg |awk 'BEGIN{OFS="\t"}
/^A/{A_name=$0}
/^B  /{B_name=$0}
/^D/{D_name=$0;print A_name,B_name,D_name}' >temp1

cat temp1|sed 's/A[0-9][0-9][0-9][0-9][0-9] //g'|sed 's/B  [0-9][0-9][0-9][0-9][0-9] //g'|sed 's/D      //g'|awk '
{FS=OFS="\t"}
{split($3,x," ");print $1,$2,x[1]}' > kegg_class.txt



cat ${group}.KEGG_Annotation.txt|cut -f 1,2,3|sort -k1,1 -u > temp

#count
awk '
{FS=OFS="\t"}
NR==FNR{a[$3]=$1"\t"$2;next}
{print $1,$2,$3,a[$2]}' kegg_class.txt temp > ${group}.KEGG_Annotation_Class.txt

cat ${group}.KEGG_Annotation_Class.txt|awk '
BEGIN{FS=OFS="\t"}
{p=$4"\t"$5;a[p]++}
END{for (i in a){print i,a[i]}}' >${group}.kegg_count.out


#表头
group=G-24vsCk-24 
#group=10-24vsCk-24
sed "1 i Gene_ID\tKEGG\tDescription\tA_Class\tB_Class" ${group}.KEGG_Annotation_Class.txt > temp
mv temp ${group}.KEGG_Annotation_Class.txt


###
#G-24 vs 10-24
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ../../R_input/CountMatrix4DESeq.csv ../../R_input/colData.csv 10-24 ./ FALSE

