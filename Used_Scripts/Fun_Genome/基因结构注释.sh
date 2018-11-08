#data prep

# Cardamine hirsutat基因组数据
mkdir chi_annotation && cd chi_annotation
wget http://chi.mpipz.mpg.de/download/sequences/chi_v1.fa
cat chi_v1.fa | tr 'atcg' 'ATCG' > chi_unmasked.fa #??

# 注释结果
wget http://chi.mpipz.mpg.de/download/annotations/carhr38.gff

# Cardamine hirsutat转录组数据
mkdir rna-seq && cd rna-seq
wget -4 -q -A '*.fastq.gz' -np -nd -r 2 http://chi.mpipz.mpg.de/download/fruit_rnaseq/cardamine_hirsuta/ &
wget -4 -q -A '*.fastq.gz' -np -nd -r 2 http://chi.mpipz.mpg.de/download/leaf_rnaseq/cardamine_hirsuta/ &


##########------------1.重复序列屏蔽-------------###########
mkdir 00-RepeatMask
RepeatMasker -e ncbi -species arabidopsis -pa 5 -gff -dir 00-RepeatMask/ chi_annotation/chi_unmasked.fa
# -e ncbi
# -species 选择物种 用/data1/software/miniconda2/share/RepeatMasker/util/queryRepeatDatabase.pl -tree 了解
# -lib 增加额外数据库,
# -pa 并行计算
# -gff 输出gff注释
# -dir 输出路径
# annotation with the library produced by RepeatMode

# de novo predict--RepeatModel
~/opt/biosoft/RepeatModeler-open-1.0.11/BuildDatabase -name test -engine ncbi output.fa
~/opt/biosoft/RepeatModeler-open-1.0.11/RepeatModeler -database test


##########------------2.从头预测-------------#############
####augustus：依赖已有的模型
conda create -n annotation augustus=3.3
source activate annotation

seqkit faidx TAIR10.fa Chr1:1-8000 > test.fa
augustus --speices=arabidopsis test.fa > test.gff
#已经被训练的物种信息可以用augustus --species=help查看。

mkdir 01-augustsus && cd 01-augustsus
ln ../00-RepeatMask/chi_unmasked.fa.masked genome.fa
seqkit split genome.fa 
#结果文件在genome.fa.split
find genome.fa.split/ -type f -name "*.fa" | parallel -j 30 augustus --species=arabidopsis --gff3=on >> temp.gff 
#并行处理
join_aug_pred.pl < temp.gff  | grep -v '^#' > temp.joined.gff
bedtools sort -i temp.joined.gff > augustsus.gff



#genemark-ET：唯一一款支持无监督训练模型，之后再识别真核基因组蛋白编码区的工具。
gmes_petap.pl --ES --sequence genome.fa --cores 50
#genemark.gtf，是标准的GTF格式，可以使用Sequence Ontology Project提供的gtf2gff3.pl进行转换
wget http://genes.mit.edu/burgelab/miso/scripts/gtf2gff3.pl
chmod 755 gtf2gff3.pl
gtf2gff3.pl genemark.gtf | bedtools sort -i - > genemark.gff

#BRAKER：GeneMark-ET根据RNA-seq无监督训练模型寻找基因，然后用AUGUSTUS进行模型训练，最后完成基因预测
# 项目根目录
mkdir index
hisat2-build 02-augustus/genome.fa index/chi_masked
hisat2 -p 15 -x index/chi_masked -1 rna-seq/leaf_ox_r1_1.fastq.gz -2 rna-seq/leaf_ox_r1_2.fastq.gz | samtools sort -@ 10 > 03-braker/leaf_ox_r1.bam &
hisat2 -p 15 -x index/chi_masked -1 rna-seq/ox_flower9_rep1_1.fastq.gz -2 rna-seq/ox_flower9_rep1_2.fastq.gz | samtools sort -@ 10 > 03-braker/ox_flower9.bam &
hisat2 -p 15 -x index/chi_masked -1 rna-seq/ox_flower16_rep1_1.fastq.gz -2 rna-seq/ox_flower16_rep1_2.fastq.gz | samtools sort -@ 10 > 03-braker/ox_flower16.bam &


#augustus是annotation环境中 
#以未屏蔽重复序列的参考序列和BAM文件作为输入
#braker.pl --gff3 --cores 10 --species=carhr --genome=chi_unmasked.fa --bam=02-barker/leaf_ox_r1.bam,02-barker/ox_flower16.bam,02-barker/ox_flower9.bam
source activate annotation
braker.pl --gff3 --cores 16 --species=carhr --genome=chi_annotation/chi_unmasked.fa --bam=03-braker/leaf_ox_r1.bam
# --gff3: 输出GFF3格式
# --genome: 基因组序列 不用masked
# --bam: 比对后的BAM文件，允许多个
# --cores: 处理核心数
# --species:随便定义 作为输出文件夹名字

#输出文件夹叫braker 


##########------------3.同源预测-------------###########

#选择一段序列作为同源序列
seqkit faidx chi_unmasked.fa Chr1:1-5000 > chr1_5k.fa
#genome threader
gth -genomic chi_annotation/chr1_5k.fa -protein chi_annotation/cer.fa -intermediate -gff3out  > gth.gff3
# 其中cer.fa就是AT1G02205.2的氨基酸序列

#不同物种的同源注释结果
#run seperately
gth -species arabidopsis -translationtable 1 -gff3 -intermediate -protein ~/db/protein_db/Athaliana_167_TAIR10.protein.fa.gz -genomic chi_unmasked.fa -o 03-genomethreader/Athaliana.gff3 
gth -species arabidopsis -translationtable 1 -gff3 -intermediate -protein ~/db/protein_db/BrapaFPsc_277_v1.3.protein.fa.gz -genomic chi_unmasked.fa -o 03-genomethreader/Brapa.gff3 
gth -species arabidopsis -translationtable 1 -gff3 -intermediate -protein ~/db/protein_db/Osativa_323_v7.0.protein.fa.gz -genomic chi_unmasked.fa -o 03-genomethreader/Osativa.gff3 


##########------------4.转录组注释-------------###########

cd rna-seq
Trinity --seqType fq --CPU 15 --max_memory 32G --left leaf_ox_r1_1.fastq.gz,ox_flower16_rep1_1.fastq.gz,ox_flower9_rep1_1.fastq.gz --right leaf_ox_r1_2.fastq.gz,ox_flower16_rep1_2.fastq.gz,ox_flower9_rep1_2.fastq.gz 


#PASA配置

cp /data1/software/PASApipeline/pasa_conf/pasa.alignAssembly.Template.txt alignAssembly.config

# 修改如下内容
DATABASE=database.sqlite
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=80
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=9

/data1/software/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g ../chi_annotation/chi_unmasked.fa -t ../rna-seq/trinity_out_dir/Trinity.fasta --ALIGNERS blat,gmap

#中断
/data1/softwafre/PASApipeline/scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -s 2 -R -g ../chi_annotation/chi_unmasked.fa -t ../rna-seq/trinity_out_dir/Trinity.fasta --ALIGNERS blat,gmap


##########------------5.EVM整合结果-------------###########

mkdir 05-EVM && cd 05-EVM

#1.创建权重文件 
vim weights.txt

ABINITIO_PREDICTION	augustus	4
TRANSCRIPT	assembler-database.sqlite	10
#2.分割原始数据
ln ../03-braker/braker/carhr/augustus.hints.gff3 gene_predictions.gff3
ln ../04-PASA/database.sqlite.pasa_assemblies.gff3 transcript_alignments.gff3
partition_EVM_inputs.pl --genome ../chi_annotation/chi_unmasked.fa --gene_predictions gene_predictions.gff3 --transcript_alignments transcript_alignments.gff3 --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
#3.并行执行
write_EVM_commands.pl --genome ../chi_annotation/chi_unmasked.fa --weights `pwd`/weights.txt \
--gene_predictions gene_predictions.gff3 \
--transcript_alignments transcript_alignments.gff3 \
--output_file_name evm.out  --partitions partitions_list.out >  commands.list
parallel --jobs 10 < commands.list
#4.合并并行结果
recombine_EVM_partial_outputs.pl --partitions partitions_list.out --output_file_name evm.out
#5.结果转换成gff3
convert_EVM_outputs_to_GFF3.pl  --partitions partitions_list.out --output evm.out  --genome ../chi_annotation/chi_unmasked.fa
find . -regex ".*evm.out.gff3" -exec cat {} \; | bedtools sort -i - > EVM.all.gff


