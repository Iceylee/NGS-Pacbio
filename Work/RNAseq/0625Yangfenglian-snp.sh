# samtools vcf  每个文件
# variants calling by samtools
samtools mpileup -t AD,ADF,ADR,DP,SP -ugf cDNA.unique.fa WR180589S_clean_sorted.bam | bcftools call -vm  > WR180589S.samtools.vcf

##两两分析

# GATK concordance
# 选取GATK和samtools一致的结果 取交集
java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R  /data2/ClientData/2018_06/YangFengLian/Genome/cDNA.unique.fa --variant WR180593S.samtools.vcf --concordance WR180594S.samtools.vcf -o FB_tmp.samtools.vcf

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R  /data2/ClientData/2018_06/YangFengLian/Genome/cDNA.unique.fa --variant FB_tmp.samtools.vcf --concordance WR180595S.samtools.vcf -o FD_tmp.samtools.vcf


for file in `ls *_clean_sorted.bam`
do
	sample=${file/_clean_sorted.bam/}
	echo "samtools mpileup -t AD,ADF,ADR,DP,SP -ugf /data2/ClientData/2018_06/YangFengLian/Genome/cDNA.unique.fa $file | bcftools call -vm  > ../5.SNP_Indel/$sample.samtools.vcf"
done > snp.cmd



java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R  /data2/ClientData/2018_06/YangFengLian/Genome/cDNA.unique.fa --variant WR180608S.samtools.vcf --concordance WR180609S.samtools.vcf -o MC_tmp.samtools.vcf
java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R  /data2/ClientData/2018_06/YangFengLian/Genome/cDNA.unique.fa --variant MC_tmp.samtools.vcf --concordance WR180610S.samtools.vcf -o MC.samtools.vcf


#dict file
java -jar picard.jar CreateSequenceDictionary R=./resources/Homo_sapiens_assembly18.fasta O=./resources/Homo_sapiens_assembly18.dict
samtools faidx ./resources/Homo_sapiens_assembly18.fasta


#SNP注释（执行目录 5.SNP_Indel）：
for sample in MB MD MC
do
	echo "nohup java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s Annotation/${sample}.SnpEff.html -c snpEff.config  -v -ud 500 Sitotroga ${sample}.samtools.vcf > Annotation/${sample}.SnpEff.csv &"
done


#参数 Sitotroga 为构建的基因组
cd /data1/software/snpEff/data
mkdir Sitotroga #genes.gtf  sequences.fa  snpEffectPredictor.bin

vim snpEff.config
#add following in """
 """
# Sitotroga genome, version cDNAv1
Sitotroga.genome : Sitotroga
"""

#在snpEff目录
java -jar snpEff.jar build -gff3 -v Sitotroga














