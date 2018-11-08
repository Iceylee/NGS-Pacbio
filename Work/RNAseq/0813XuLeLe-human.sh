'''
1.NC实验组 SH1，SH2对照组
NCvsSH1   NCvsSH2
差异分析要跑两次 之后重跑lncRNA分析
2.无生物学重复 Gfold 单独做MA图
3.SNP分析中




'''



bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_08/xulele  /data2/ClientData/2018_08/XuLeLe

python TransSeq_WithGenome.py TransSeq_WithGenome.conf 1>0814.log 2>&1

python TransSeq_WithGenome2.py TransSeq_WithGenome.conf 1>0816.log 2>&1
#

samtools view -u WR181134S_clean.sam | samtools sort -@ 15 - > WR181134S_clean_sorted.bam

gtf=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf
htseq-count -f sam WR181134S_clean.sam $gtf -q > WR181134S_clean_CountNum.txt

###补 一组
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R /data2/ClientData/2018_08/XuLeLe/R_input/CountMatrix4DESeq.csv /data2/ClientData/2018_08/XuLeLe/R_input/colData.csv SH2 /data2/ClientData/2018_08/XuLeLe/ FALSE

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa SYMBOL ENTREZID /data2/ClientData/2018_08/XuLeLe/







#Gfold 差异分析：NCvsSH1 NCvsSH2
#sample2/sample1 s1为对照
#提取count数 在第三列

awk 'BEGIN{FS=",";OFS="\t"}{print $1,"1",$2,"1","1"}'  ../../R_input/CountMatrix4DESeq.csv|sed '1d' >NC.read_cnt

awk 'BEGIN{FS=",";OFS="\t"}{print $1,"1",$3,"1","1"}'  ../../R_input/CountMatrix4DESeq.csv|sed '1d' >SH1.read_cnt

awk 'BEGIN{FS=",";OFS="\t"}{print $1,"1",$4,"1","1"}'  ../../R_input/CountMatrix4DESeq.csv|sed '1d' >SH2.read_cnt

## 得到的Gfold比实际要小（负）or大（正）的概率为sc （sc改变，Gfold也会随着改变）
##sc 增大，Gfold的绝对值是变大
gfold diff -s1 SH1 -s2 NC -suf .read_cnt -sc 0.05 -o NCvsSH1.diff

cat NCvsSH1.diff|grep -v "#"|awk '$3>=1||$3<=-1{print $0}'>NCvsSH1.sig.diff

##
gfold diff -s1 SH2 -s2 NC -suf .read_cnt  -sc 0.05 -o NCvsSH2.diff

cat NCvsSH2.diff|grep -v "#"|awk '$3>=1||$3<=-1{print $0}'>NCvsSH2.sig.diff


##
exp=NC
base=SH1
cat ${exp}vs${base}.diff|grep -v "#"|cut -f 1,3 > temp

#第2，3列求平均
cat CountMatrix4DESeq.csv|awk 'BEGIN{FS=",";OFS="\t"}{$1=$1;print $1,($2+$3)/2}' > NCvsSH1.countMean
cat CountMatrix4DESeq.csv|awk 'BEGIN{FS=",";OFS="\t"}{$1=$1;print $1,($2+$4)/2}' > NCvsSH2.countMean


#根据geneID提取count平均值
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2;next}
{if($2 <=-1 || $2 >=1)sig="TRUE";
if($2>-1 && $2<1)sig="FALSE";
	print $1,a[$1],$2,sig}' ${exp}vs${base}.countMean temp |sed "1 i geneID\tcountMean\tlog2FoldChange\tSignificant"> ${exp}vs${base}_all_genes_exprData.txt

#sig gene
cat ${exp}vs${base}_all_genes_exprData.txt|grep TRUE|sed "1 i geneID\tcountMean\tlog2FoldChange\tSignificant" > ${exp}vs${base}_sig_genes_exprData.txt









#Result
mkdir 0.QC 1.Mapping 2.GenesExpress 3.DiffExprGene 4.GO_KEGG_Enrichment  5.SNP_Indel  6.lncRNA_Expr 7.Correlation_analysis

ln -s ../../5.SNP_Indel/*6.gatk.vcf ./
ln -s ../../5.SNP_Indel/*SNPs.vcf ./
ln -s ../../5.SNP_Indel/*SNPs.PASS.vcf ./
ln -s ../../5.SNP_Indel/*INDELs.PASS.vcf ./
ln -s ../../5.SNP_Indel/*INDELs.vcf ./
ln -s ../../5.SNP_Indel/*AllInfo.txt ./
ln -s ../../5.SNP_Indel/*SnpEff.html ./
ln -s ../../5.SNP_Indel/*SnpEff.genes.txt ./











