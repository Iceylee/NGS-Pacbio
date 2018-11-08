#0.1 NCBI下载的sra格式
wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/SRR177/SRR1770413/SRR1770413.sra
#sratoolkit 得到read1和read2的fastq
fastq-dump --split-files SRR1770413.sra
#bgzip（不推荐gzip） 压缩
bgzip -f SRR1770413_1.fastq
bgzip -f SRR1770413_2.fastq
#gzip SRR1770413_1.fastq > SRR1770413_1.fastq.gz


#0.2 fasta
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna E.coli_K12_MG1655.fa
samtools faidx E.coli_K12_MG1655.fa
#获取一段序列
samtools faidx E.coli_K12_MG1655.fa NC_000913.3:1000000-1000200



#1.比对
#建立比对索引
bwa index E.coli_K12_MG1655.fa

#1.1 比对
time bwa mem -t 4 -R '@RG\tID:foo\tPL:illumina\tSM:E.coli_K12' input/E.coli/fasta/E.coli_K12_MG1655.fa input/E.coli/fastq/SRR1770413_1.fastq.gz input/E.coli/fastq/SRR1770413_2.fastq.gz | samtools view -Sb - > output/E.coli/E_coli_K12.bam && echo "** bwa mapping done **"
#-R 设置Read Group信息，它是read数据的组别标识，并且其中的ID，PL和SM信息在正式的项目中是不能缺少的(如果样本包含多个测序文库的话，LB信息也不要省略)，另外由于考虑到与GATK的兼容关系，PL（测序平台）信息不能随意指定，必须是：ILLUMINA，SLX，SOLEXA，SOLID，454，LS454，COMPLETE，PACBIO，IONTORRENT，CAPILLARY，HELICOS或UNKNOWN这12个中的一个。


#1.2 按参考序列位置的从小到大排序
time samtools sort -@ 4 -m 4G -O bam -o output/E.coli/E_coli_K12.sorted.bam output/E.coli/E_coli_K12.bam && echo "** BAM sort done"
rm -f output/E.coli/E_coli_K12.bam

#1.3 使用GATK标记出排完序的数据中的PCR重复序列。
time gatk MarkDuplicates -I output/E.coli/E_coli_K12.sorted.bam -O output/E.coli/E_coli_K12.sorted.markdup.bam -M output/E.coli/E_coli_K12.sorted.markdup_metrics.txt && echo "** markdup done **"


#1.4 删除不必要文件(可选)
rm -f output/E.coli/E_coli_K12.bam
rm -f output/E.coli/E_coli_K12.sorted.bam

#1.5 创建比对索引文件
time samtools index output/E.coli/E_coli_K12.sorted.markdup.bam && echo "** index done **"


# 2.变异检测
#在input的fasta文件夹下
#先为E.coliK12的参考序列生成一个.dict文件，这可以通过调用CreateSequenceDictonary模块来完成(这是原来picard的功能)。
#.dict文件的名字前缀需要和fasta的一样，并跟它在同一个路径下，这样GATK才能够找到。
gatk CreateSequenceDictionary -R E.coli_K12_MG1655.fa -O E.coli_K12_MG1655.dict && echo "** dict done **"

#先为每个样本生成一个GVCF，然后再用GenotypeGVCFs对这些GVCF进行joint calling
#1 生成中间文件gvcf
time gatk HaplotypeCaller \
 -R input/E.coli/fasta/E.coli_K12_MG1655.fa \
 --emit-ref-confidence GVCF \
 -I output/E.coli/E_coli_K12.sorted.markdup.bam \
 -O output/E.coli/E_coli_K12.g.vcf && echo "** gvcf done **"

#2 通过gvcf检测变异
time gatk GenotypeGVCFs \
 -R input/E.coli/fasta/E.coli_K12_MG1655.fa \
 -V output/E.coli/E_coli_K12.g.vcf \
 -O output/E.coli/E_coli_K12.vcf && echo "** vcf done **"

#获得了E.coli K12这个样本初步的变异结果——E_coli_K12.vcf。
#之所以非要分成两个步骤，是因为我想借此告诉大家，变异检测不是一个样本的事情，有越多的同类样本放在一起joint calling结果将会越准确，而如果样本足够多的话，在低测序深度的情况下也同样可以获得完整并且准确的结果，而这样的分步方式是应对多样本的好方法。

# 压缩 
time /Tools/common/bin/bgzip -f output/E.coli/E_coli_K12.vcf
# 构建tabix索引
time /Tools/common/bin/tabix -p vcf output/E.coli/E_coli_K12.vcf.gz



#重比对和BQSR
#VQSR 变异质控 过滤变异结果

#硬过滤
# 使用SelectVariants，选出SNP
time gatk SelectVariants \
    -select-type SNP \
    -V ../output/E.coli/E_coli_K12.vcf.gz \
    -O ../output/E.coli/E_coli_K12.snp.vcf.gz

# 为SNP作硬过滤
time gatk VariantFiltration \
    -V ../output/E.coli/E_coli_K12.snp.vcf.gz \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ../output/E.coli/E_coli_K12.snp.filter.vcf.gz

# 使用SelectVariants，选出Indel
time gatk SelectVariants \
    -select-type INDEL \
    -V ../output/E.coli/E_coli_K12.vcf.gz \
    -O ../output/E.coli/E_coli_K12.indel.vcf.gz

# 为Indel作过滤
time gatk VariantFiltration \
    -V ../output/E.coli/E_coli_K12.indel.vcf.gz \
    --filter-expression "QD < 2.0 || FS > 200.0 || SOR > 10.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "PASS" \
    -O ../output/E.coli/E_coli_K12.indel.filter.vcf.gz

# 重新合并过滤后的SNP和Indel
time gatk/4.0.1.2/gatk MergeVcfs \
    -I ../output/E.coli/E_coli_K12.snp.filter.vcf.gz \
    -I ../output/E.coli/E_coli_K12.indel.filter.vcf.gz \
    -O ../output/E.coli/E_coli_K12.filter.vcf.gz

# 删除无用中间文件
rm -f ../output/E.coli/E_coli_K12.snp.vcf.gz* ../output/E.coli/E_coli_K12.snp.filter.vcf.gz* ../output/E.coli/E_coli_K12.indel.vcf.gz* ../output/E.coli/E_coli_K12.indel.filter.vcf.gz*







