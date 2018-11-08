################################################
##0.install
################################################

docker pull wwliao/cnvnator
conda install cnvkit

BL实验 LM对照







################################################
##1.QC & trim
################################################

bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_06/SongJianCheng/Macrogen/LM /data2/ClientData/2018_06/SongJianCheng/clean_data

bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_06/SongJianCheng/Macrogen/BL /data2/ClientData/2018_06/SongJianCheng/clean_data







################################################
##2.Mapping
################################################

##################
##2.1 bwa
##################
#bwa index  -p genome 即前缀为genome
cd /data1/GenomicDatabases/Maize
bwa index GCF_000005005.2_B73_RefGen_v4_genomic

#bwa map
for smp in BL LM
do
	ref=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
	fq1=../0.QC/clean_data/${smp}_clean_R1.fastq.gz
	fq2=../0.QC/clean_data/${smp}_clean_R2.fastq.gz
	echo "bwa mem -t 20 $ref $fq1 $fq2 > ${smp}.sam"
done >cmd


samtools view -u LM.sam | samtools sort -@ 15 - > LM_sorted.bam 
samtools view -u BL.sam | samtools sort -@ 15 - > BL_sorted.bam


samtools depth LM_sorted.bam > LM_depth.txt
samtools depth BL_sorted.bam > BL_depth.txt


bedtools genomecov -ibam LM_sorted.bam -d >LM.dep_base
bedtools genomecov -ibam BL_sorted.bam -d >BL.dep_base

################################################
##4.SV
################################################
breakdancer-max


breakdancer_DIR=/data1/software/miniconda2/pkgs/breakdancer-1.4.5-2/bin
#Coefficient of variation *** in library *** is larger than the cutoff 1, poor quality data, excluding from further analysis.  change cutoff  -v 3 
perl $breakdancer_DIR/bam2cfg.pl ../1.Mapping/BL_sorted.bam > BL.cfg #1min
perl $breakdancer_DIR/bam2cfg.pl ../1.Mapping/LM_sorted.bam > LM.cfg
#只保留最后两行


breakdancer-max BL.cfg > BL_SV.txt #10am
breakdancer-max LM.cfg > LM_SV.txt 

#-q 35 
# -m 100 -x 1000000 -s 30 -d 5

#DEL (缺失 )1, INS (插入 )，INV (倒位 )2， ITX (染色体内易位 )3， CTX (染色体间易位)4


for smp in BL LM
do
      #去掉第11列 （每个map file中分别多少个）
	cat ${smp}_SV.txt|awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$12}' > temp

	  #filter : 1)reads >=5 $10 2)100<SV_Size<1000000 $8 3)quality >=35(默认参数)
	cat temp |awk '
	BEGIN{FS=OFS="\t"}
	$0~/#/{print $0}
	$0!~/#/ &&$10>=5 && $8>100 && $8<1000000{print $0}' > ${smp}_SV.tsv

done

#统计
cat BL_SV.tsv|awk '
{a[$7]+=1}END{for (i in a){print i,a[i]}}' > BL.stats
cat LM_SV.tsv|awk '
{a[$7]+=1}END{for (i in a){print i,a[i]}}' > LM.stats


################################################
##5.CNV
################################################

conda install pyfaidx
cd /data1/GenomicDatabases/Maize/split
faidx -x ../GCF_000005005.2_B73_RefGen_v4_genomic.fa

#cpu 可用到30-40cpu 需要限制。
sudo service docker start
docker run \
  --user root \
  -it \
  -v /data2/ClientData/2018_06/SongJianCheng/1.Mapping/:/bam \
  -v /data1/GenomicDatabases/Maize/:/ref \
  -v /data2/ClientData/2018_06/SongJianCheng/5.CNV:/data \
  --rm \
  #--cpus=20 \
  wwliao/cnvnator \
  /bin/bash



#docker exec -it ed0244b9d514 /bin/bash
#ctrl + P + Q

ref=/ref/GCF_000005005.2_B73_RefGen_v4_genomic.fa

for smp in BL LM
do
	cnvnator -genome $ref -root $smp.root -tree /bam/${smp}_sorted.bam -unique
	cnvnator -genome $ref -root $smp.root -d /ref/split/ -his 100 #d路径为参考基因fa按照染色体名字分隔开的 chr1.fa chr2.fa 所在目录
	cnvnator -genome $ref -root $smp.root -stat 100
	cnvnator -genome $ref -root $smp.root -partition 100
	cnvnator -genome $ref -root $smp.root -call 100 > $smp.cnvout.txt
done #2h

#filter:1)q0>=0.5 2)e-val1>=0.01
for smp in BL LM
do
	cat $smp.cnvout.txt |awk '
		BEGIN{FS=OFS="\t"}
		$5<0.01 && $9<0.5 {print $0}'|sed '1 i CNV_type\tcoordinates\tCNV_size\tnormalized_RD\te-val1\te-val2\te-val3\te-val4\tq0' > ${smp}_cnv.tsv
done







##cnvkit
cnvkit.py batch -n *Normal.bam --output-reference new_reference.cnn \
    -t my_targets.bed -a my_antitargets.bed --male-reference \
    -f hg19.fasta -g data/access-5kb-mappable.hg19.bed


#全基因组
cnvkit.py batch 700_bwa.sam.bam --annotate ucsc-human-refflat2.txt --normal 699_bwa.sam.bam  --method wgs -f hg38.fasta  --output-reference my_flat_reference.cnn -d 699vs700

#trial
ref=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
cnvkit.py batch BL_sorted.bam --normal LM_sorted.bam  --method wgs -f $ref  --output-reference my_flat_reference.cnn -d BLvsLM
#BL 实验

################################################
##6.snp
################################################

for smp in BL LM
do
	ref=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
	echo "nohup samtools mpileup -t AD,ADF,ADR,DP,SP -ugf $ref ${smp}_sorted.bam | bcftools call -vm  > ../3.SNP_Indel/$smp.samtools.vcf &"
done > snp.cmd

 #### Variant filtering ####

variant_name=BL.samtools.vcf
filter_name=BL.fit1.vcf
filter_log_name=BL.flt1.log
GenomeFa=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
GenomeDict=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.dict

java -jar /data1/software/picard/picard.jar CreateSequenceDictionary R=$GenomeFa O=$GenomeDict #1min

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T VariantFiltration  -R $GenomeFa -V $variant_name -window 35 -cluster 3 -filterName FSS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o $filter_name 2>$filter_log_name #3h


 #### SNP/INDEL ####
SNP_name=BL.SNPs.vcf
INDEL_name=BL.INDELs.vcf
filter_name=BL.fit1.vcf
GenomeFa=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType SNP -o $SNP_name # 1h

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType INDEL -o $INDEL_name #1h








 #### Annotation ####
#参数 Maize 为构建的基因组
cd /data1/software/snpEff/data
mkdir Maize #genes.gtf  sequences.fa  snpEffectPredictor.bin
cp /data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa Maize/sequences.fa
cp /data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.gtf Maize/genes.gtf

vim snpEff.config
#add following in """
 """
# Maize, version GCF_000005005.2_B73_RefGen_v4_genomic
Maize.genome : Maize
"""

#在snpEff目录
java -jar snpEff.jar build -gff3 -v Maize #gff
java -jar snpEff.jar build -gtf22 -v Maize #gtf


filter_name=BL.fit1.vcf
Annota_Html=BL.SnpEff.html
Annota_vcf=BL.SnpEff.vcf
snpEff_SpeciesName=Maize
java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s $Annota_Html -c /data1/software/snpEff/snpEff.config -v -ud 500 $snpEff_SpeciesName $filter_name > $Annota_vcf


####plot
ref=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
seqkit fx2tab -l -g -n -i -H $ref |cut -f 1,4 > contig_len.txt

#name   length
NC_024459.2     307041717
NC_024460.2     244442276
NC_024461.2     235667834
NC_024462.2     246994605
NC_024463.2     223902240
NC_024464.2     174033170
NC_024465.2     182381542
NC_024466.2     181122637
NC_024467.2     159769782
NC_024468.2     150982314


#取上面10chr
cat BL_depth.txt|awk '
BEGIN{FS=OFS="\t"}
$1=="NC_024459.2"||$1=="NC_024460.2"||$1=="NC_024461.2"||$1=="NC_024462.2"||$1=="NC_024463.2"||$1=="NC_024464.2"||$1=="NC_024465.2"||$1=="NC_024466.2"||$1=="NC_024467.2"||$1=="NC_024468.2" {print $0}' > BL_depth_filter.txt

#100kb window
cat BL.head.txt|awk '
BEGIN{FS=OFS="\t"}
{loc=int($2/100000);a[loc]=a[loc]+$3}
END{for (i in a) {print i,a[i]}}' 

cat BL.head.txt|awk '
BEGIN{FS=OFS="\t"}
$2%100000!=0{sum=sum+$3}
$2%100000==0{sum=sum+$3;print $1,$2,sum}' |head

cat LM_depth_filter.txt|awk '
BEGIN{FS=OFS="\t"}
$2%100000==1{sum=0}
$2%100000!=0{sum=sum+$3}
$2%100000==0{sum=sum+$3;print $1,$2,sum}' > LM_depth_filter_win.txt



for chr in NC_024459.2 NC_024460.2 NC_024461.2 NC_024462.2 NC_024463.2 
do
	cat BL_depth_filter.txt|awk -v chr="$chr" '
		BEGIN{FS=OFS="\t"}
		$1==chr {print $0}' > BL_${chr}.txt
done

for chr in NC_024464.2 NC_024465.2 NC_024466.2 NC_024467.2 NC_024468.2
do
	cat BL_depth_filter.txt|awk -v chr="$chr" '
		BEGIN{FS=OFS="\t"}
		$1==chr {print $0}' > BL_${chr}.txt
done

for chr in NC_024464.2 NC_024465.2 NC_024466.2 NC_024467.2 NC_024468.2
do
	for smp in LM
	do
		cat ${smp}_${chr}.txt|awk '
			BEGIN{FS=OFS="\t"}
			{loc=int($2/100000);a[loc]=a[loc]+$3}
			END{for (i in a) {print i,a[i]}}' > ${smp}_${chr}.win.txt
	done
done

nohup wc -l BL_depth.txt > BL.wc.txt &
#coverage:depth.txt行数即测到的总碱基数/参考基因组大小
BL:1635254242/2135083061 0.76589
LM:1808734238/2135083061 0.84715
ref:2135083061

#depth：depth.txt第三列总和即总count数/测到的总碱基数
BL：64494083992/1635254242 39.44
LM：62704846150/1808734238 34.67



nohup cat BL_depth.txt|awk '{sum=sum+$3}END{print sum}' >BL.sum &
nohup cat LM_depth.txt|awk '{sum=sum+$3}END{print sum}' >LM.sum &



################################################
##7.kegg&go
################################################

grep -v "#"  BL.SnpEff.genes.txt|cut -f 1|sort -k1,1 -u > BL.gene.list
grep -v "#"  LM.SnpEff.genes.txt|cut -f 1|sort -k1,1 -u > LM.gene.list





#打包

tar cvzf - ./picture | split -d -b 1000m - picture
tar cvzf - ./sam 
cat BL* >BL.tar.gz











################################################
##3.Fusion Gene
################################################
#db download
nohup wget https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/GRCh38_v27_CTAT_lib_Feb092018.plug-n-play.tar.gz &

#install
docker pull trinityctat/ctatfusion

# setting some targets
#VERSION=1.4.0
CTAT_GENOME_LIB="GRCh38_v27_CTAT_lib_Feb092018"
SMP=BL

# and now running STAR-Fusion & FusionInspector 'inspect' & Trinity de-novo reconstruction via Docker:
cd 2.Fusion_Gene

docker run -v `pwd`/../:/data \
    -v /data1/CTAT_Genome/GRCh38_v27_CTAT_lib_Feb092018/:/database \
    --rm --cpus=15 trinityctat/ctatfusion \
    /usr/local/src/STAR-Fusion/STAR-Fusion \
    --left_fq /data/0.QC/clean_data/BL_clean_R1.fastq.gz \
    --right_fq /data/0.QC/clean_data/BL_clean_R2.fastq.gz \
    --genome_lib_dir /database/ctat_genome_lib_build_dir \
    -O /data/BL_StarFusionOut \
    --FusionInspector validate \
    --examine_coding_effect \
    --denovo_reconstruct

#
star --runThreadN  5 --genomeDir $hg19_star_index --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  \
--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 \
--chimJunctionOverhangMin 12  --alignSJDBoverhangMin 10  --alignMatesGapMax 100000 \
--alignIntronMax 100000 --chimSegmentReadGapMax parameter 3  --alignSJstitchMismatchNmax 5 -1 5 5 \
--readFilesIn  $fq1 $fq2 --outFileNamePrefix  ${sample}_star 

STAR-Fusion --genome_lib_dir /data1/CTAT_Genome/GRCh38_v27_CTAT_lib_Feb092018/ctat_genome_lib_build_dir/   -J std.Chimeric.out.junction  --output_dir star_fusion_outdir



##################
##2.1 STAR
##################
#STAR index
cd /data1/GenomicDatabases/Maize

STAR --runThreadN 15 --runMode genomeGenerate --genomeDir star_index/ --genomeFastaFiles GCF_000005005.2_B73_RefGen_v4_genomic.fa #10min

#STAR mapping


for smp in BL LM
do
	ref=/data1/GenomicDatabases/Maize/star_index/
	fq1=../0.QC/clean_data/${smp}_clean_R1.fastq.gz
	fq2=../0.QC/clean_data/${smp}_clean_R2.fastq.gz
	echo "STAR --runThreadN  15 --genomeDir $ref --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  \
--twopassMode Basic --outReadsUnmapped None --chimSegmentMin 12 \
--chimJunctionOverhangMin 12  --alignSJDBoverhangMin 10  --alignMatesGapMax 100000 \
--alignIntronMax 100000 --chimSegmentReadGapMax parameter 3  --alignSJstitchMismatchNmax 5 -1 5 5 \
--readFilesIn  $fq1 $fq2 --outFileNamePrefix  ${smp}_star --limitSjdbInsertNsj 2000000"
done >cmd


