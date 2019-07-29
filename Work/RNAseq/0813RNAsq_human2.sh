'''
1.两组6个样本 人源
2.差异基因较少
3.前50聚类热图去掉了第一个基因（异常表达）

'''

#0.QC

#Preparation
sudo chmod g+w -R XieHuaPing/

#1.参考基因组
# mkdir Genome
# cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.fa Genome/
# cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.gtf Genome/

#2.pipeline
cp /data1/script/Pipeline/RNA-Seq/* ./

#3.files needed
mkdir R_input
mv colData.csv R_input/colData.csv.temp
mv sample_list.txt sample_list.txt.temp

LZY-Con-1	WR181074S
LZY-Con-2	WR181075S
LZY-Con-3	WR181076S
LZY-MA-1	WR181077S
LZY-MA-2	WR181078S
LZY-MA-3	WR181079S
exp:LZY-MA
对照：LZY-Con










#conf

%s/\/data3\/ClientData\/2018_07\//\/data2\/ClientData\/2018_08\//g
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/0Get_OrgDb_KEGG_Name.R “Danio rerio”

#[1] "AH57981"
#   kegg_code scientific_name common_name
# 71       dre     Danio rerio   zebrafish

#snpeff
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.fa sequences.fa
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.gtf genes.gtf

java -jar snpEff.jar build -gtf22 -v Danio_rerio #gtf

#修改名字
for smp in WR181061S WR181062S
do
	mv ${smp}-1_clean_R1.fastq.gz ${smp}_1_clean_R1.fastq.gz
	mv ${smp}-1_clean_R2.fastq.gz ${smp}_1_clean_R2.fastq.gz
done



samtools view -u WR181077S_clean.sam | samtools sort -@ 15 - > WR181077S_clean_sorted.bam

gtf=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf
htseq-count -f sam WR181077S_clean.sam $gtf -q > WR181077S_clean_CountNum.txt


Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa SYMBOL ENTREZID /data2/ClientData/2018_08/XieHuaPing/

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa SYMBOL ENTREZID /data2/ClientData/2018_08/XieHuaPing/





