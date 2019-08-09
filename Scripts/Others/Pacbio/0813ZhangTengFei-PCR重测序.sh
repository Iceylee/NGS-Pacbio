'''
0.原始PCR序列400多bp。基因为1.5MB左右。绵羊MSTN
1.clean_data 为我们自己过滤数据（min36） clean_data2为诺和过滤（min250） 分析用我们过滤的clean数据来做。（QC数据质量不高，我们的相对好一些）
2.bwa 和chr2进行比对。之后SNP和注释
3.序列拼接，之后cd-hit去冗余
NC_019459.2:118140420-118145410
4.GATK做snp分析结果不在目标基因范围。
改用samtools做。

'''
#0.QC

bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_08/zhengtengfei/1.rawdata  /data2/ClientData/2018_08/ZhangTengFei

#nova QC
/data1/software/FastQC/fastqc -o FastQC4Nova -t 20 clean_data2/WR180807S_clean_R1.fastq.gz clean_data2/WR180807S_clean_R2.fastq.gz

md5sum clean_data2/WR180807S_clean_R1.fastq.gz >> FastQC4Nova/clean_data.md5
md5sum clean_data2/WR180807S_clean_R2.fastq.gz >> FastQC4Nova/clean_data.md5

cd FastQC4Nova
python /data1/script/Q30_analy.py . clean

#1.alignment
 #genome : /data1/GenomicDatabases/Sheep
# gene : https://www.ncbi.nlm.nih.gov/gene/?term=443449
cat /data1/GenomicDatabases/Sheep/GCF_000298735.2_Oar_v4.0_genomic.gtf|grep MSTN > MSTN.gtf

#2.snp
snpEff_abr=MSTN_gene
python SNP_Indels.py Genome/chr2.fa ./ $snpEff_abr FALSE


#snpeff 构建基因组
#参数 Locusta_migratoria 为构建的基因组
python /data1/script/ExtraGeneID2Fa.py chr2.txt /data1/GenomicDatabases/Sheep/GCF_000298735.2_Oar_v4.0_genomic.fa > chr2.fa
chr2.gtf


##################
##2.1 bwa
##################
#bwa index  -p genome 即前缀为genome


#bwa map

ref=../../Genome/chr2.fa
fq1=../../clean_data2/WR180807S_clean_R1.fastq.gz
fq2=../../clean_data2/WR180807S_clean_R2.fastq.gz
bwa index $ref
bwa mem -t 20 $ref $fq1 $fq2 > WR180807S.sam

smp=WR180807S

samtools view -u ${smp}.sam | samtools sort -@ 15 - > ${smp}_sorted.bam 
samtools depth ${smp}_sorted.bam > ${smp}_depth.txt
#bedtools genomecov -ibam ${smp}_sorted.bam -d >${smp}.dep_base

samtools flagstat ${smp}_sorted.bam> ${smp}.stat

##################
##2.2 hisat2
##################
ref=../../Genome/NC_019459.2
fq1=../../clean_data/WR180807S_clean_R1.fastq.gz
fq2=../../clean_data/WR180807S_clean_R2.fastq.gz

hisat2 -p 15 -x $ref -1 $fq1 -2 $fq2 -S WR180807S.sam > log 2>&1




####400bp序列，PE250已测通
###1.R1和R2拼接

source activate qiime1

fq1=../../clean_data/WR180807S_clean_R1.fastq.gz
fq2=../../clean_data/WR180807S_clean_R2.fastq.gz
join_paired_ends.py -f $fq1 -r $fq2  -o WR180807S_join

###2.
seqret -sequence fastqjoin.join.fastq -outseq WR180807S_join.fa

###3.统计
seqkit fx2tab -l -g -n -i -H WR180807S_join.fa |cut -f 1,4|sort -k2nr,2 > length.txt

cat length.txt|awk 'BEGIN{FS=OFS="\t"}{loc=int($2/50);a[loc]=a[loc]+1}
   END{for (i in a)print i*50,a[i],a[i]/116058}' > stat.txt
##nova 400 以上占99%
## BYZH 400 以上占90%

cat length.txt|awk '$2>=400{print $0}' > length_400.txt

cat length_400.txt|cut -f 1 > 400.id.txt #115032

###4.仅用400以上长的序列
cat WR180807S_join.fa|seqkit grep -f 400.id.txt > 400.join.fa

#rename
cat 400.join.fa|awk '
/>/{num++;print ">Cluster"num}
!/>/{print $0}' > 400.join.rename.fa


cd-hit-est -i 400.join.rename.fa -o WR180807S_join_CDHIT.fa  -n 8 -p 1 -g 1 -M 0 -T 20 -d 0 -c 1



#统计序列
#如果以>开头，打印该列和前一列。最后 打印最后一列。
cat WR180807S_join_CDHIT.fa.clstr|awk 'BEGIN { FS=OFS="\t" }
{
    if (match($0,/^>/)!=0){
        print prev
        print $0
    }
    prev = $0
} END{print $0}'> cdhit-1.stat

cat cdhit-1.stat|awk 'BEGIN{FS=OFS="\t"}
	NR%2==0{prev=$0}
	NR%2!=0{print prev,$1+1,$2}'|sort -k3nr,3 > cdhit-1.tab


seqkit stat 400.join.rename.fa #106721
cat cdhit-1.tab|cut -f 2,3,4|awk -F '[\t,.>]' 'BEGIN{OFS="\t"}{print $4,$2,$1,$1/1067.21}'|sed '1 i ClusterID\tLength\tReads_Num\tPercentage(%)' > sequence_stat.txt

#提取序列
cat sequence_top50_stat.txt|cut -f 1 > top50.id
cat cluster.fa|seqkit grep -f top50.id > cluster_top50.fa

#plot coverage
cat WR180807S_depth.txt|awk '
			BEGIN{FS=OFS="\t"}
			{loc=int($2/100);a[loc]=a[loc]+$3}
			END{for (i in a) {print i,a[i]}}'|sort -k1nr,1 > WR180807S.win.txt

cat WR180807S_depth.txt|awk '
			BEGIN{FS=OFS="\t"}
			$2>115000000&&$2<120000000{print $0}'|sort -k2nr,1 > WR180807S_118.txt


#mean depth
 wc -l WR180807S_depth.txt #6945
cat WR180807S_depth.txt|awk '{sum=sum+$3}END{print sum/6945}'


#samtools snp
ref=../Genome/chr2.fa
bam=../1.Mapping/bwa_chr2_2/WR180807S_sorted.bam
vcf=WR180807S.samtools.vcf
smp=WR180807S

samtools mpileup -t AD,ADF,ADR,DP,SP -ugf $ref $bam | bcftools call -vm  > $vcf

vcfutils.pl varFilter -Q 20 -d 4 $vcf > ${smp}.flt1.vcf

bcftools view -v snps ${smp}.flt1.vcf > ${smp}.SNPs.vcf
bcftools view -v indels ${smp}.flt1.vcf > ${smp}.INDELs.vcf


java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s ${smp}.SnpEff.html -c /data1/software/snpEff/snpEff.config -v -ud 500 MSTN_gene ${smp}.flt1.vcf > ${smp}.SnpEff.vcf

#0822 售后 序列统计：不设过滤条件
#rename
cat WR180807S_join.fa|awk '
/>/{num++;print ">Cluster"num}
!/>/{print $0}' > WR180807S_join.rename.fa
mv WR180807S_join.rename.fa Cluster_all.fa

cd-hit-est -i Cluster_all.fa -o WR180807S_join_CDHIT.fa  -n 8 -p 1 -g 1 -M 0 -T 20 -d 0 -c 1

#统计序列
#如果以>开头，打印该列和前一列。最后 打印最后一列。
cat WR180807S_join_CDHIT.fa.clstr|awk 'BEGIN { FS=OFS="\t" }
{
    if (match($0,/^>/)!=0){
        print prev
        print $0
    }
    prev = $0
} END{print $0}'> cdhit-1.stat

cat cdhit-1.stat|awk 'BEGIN{FS=OFS="\t"}
	NR%2==0{prev=$0}
	NR%2!=0{print prev,$1+1,$2}'|sort -k3nr,3 > cdhit-1.tab

cat cdhit-1.tab|cut -f 2,3,4|awk -F '[\t,.>]' 'BEGIN{OFS="\t"}{print $4,$2,$1,$1/1067.21}'|sed '1 i ClusterID\tLength\tReads_Num\tPercentage(%)' > sequence_stat.txt


mv sequence_stat.txt sequence_stat_all.txt
