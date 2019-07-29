################################################
##1.Trinity assembly
################################################

##################
##1.1 Trinity Comstruction
##################

cd 1.Mapping
/data1/software/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --CPU 20 --max_memory 100G --samples_file sample.list --output 1.Trinity/
#11pm-5am 18h


### CD-HIT 0.9
mkdir CD-HIT;cd CD-HIT
cd-hit-est -i ../Trinity.fasta -o Trinity_CD-HIT_0.9.fa -c 0.9 -n 8 -p 1 -g 1 -M 0 -T 20 -d 0 #5min

cd-hit-est -i ../../Trinity.fasta -o Trinity_CD-HIT_0.9.fa -c 0.95 -n 10 -M 16000 -T 20 -d 0 

/data1/software/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl $fa > N50.txt 

##################
##1.2 Trinity Mapping
##################

/data1/software/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts Trinity_CD-HIT_0.9.fa  --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

#hisat2 比对
for file in `ls *_R1.fastq.gz`
do
	index=../1.Trinity/CD-HIT/index/fa
	filename=${file/_R1.fastq.gz/}
	samfile=../7.hisat2/$filename.sam
	hisat2_log=../7.hisat2/${filename}_hisat2.log
	echo "hisat2 -p 15 -x $index -1 ${filename}_R1.fastq.gz -2 ${filename}_R2.fastq.gz -S $samfile > ${hisat2_log} 2>&1"
done >cmd

hisat2-build Trinity_CD-HIT_0.9.fa index/fa

for file in WR180739S WR180740S WR180741S WR180742S
do
	echo "samtools view -u ${file}_clean.sam | samtools sort -@ 15 - > ${file}_sorted.bam"
done > sort.cmd #10am-

python /data1/script/CountBamState.py ./

for i in $(find . -name '*_hisat2.log')
do 
	echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> mapping_stat_list.txt
done

python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > Statistic_Mapping.txt

### mapping

mkdir 3.Mapping; cd 3.Mapping

for i in WR180739S WR180740S WR180741S WR180742S
do
        left=../clean_data/${i}_clean_R1.fastq.gz
        right=../clean_data/${i}_clean_R2.fastq.gz
        fa=../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa

        echo "/opt/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --thread_count 20 --transcripts $fa  --seqType fq --left $left --right $right --est_method RSEM --aln_method bowtie --trinity_mode --SS_lib_type RF --gene_trans_map $i --output_dir ${i}_rsem_outdir"

done > cmd









################################################
##2.TransDecoder CDS prediction
################################################

cd 1.Mapping

mkdir TransDecoder ; cd TransDecoder

/data1/software/TransDecoder/TransDecoder.LongOrfs  -t ../CD-HIT/Trinity_CD-HIT_0.9.fa #5min

hmmscan --cpu 8 --domtblout pfam.domtblout /data1/Pfam/Pfam-A.hmm Trinity_CD-HIT_0.9.fa.transdecoder_dir/longest_orfs.pep #9:52-6pm 

/data1/software/TransDecoder/TransDecoder.Predict -t ../CD-HIT/Trinity_CD-HIT_0.9.fa --retain_pfam_hits pfam.domtblout --cpu 20 #5min










################################################
##3.Gene Annotation 
################################################

cd 2.Annotation

##################
##3.1 nr,uniprot,kog
##################
pep=../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa
out=nr.out
dbdir=/data1/DiamondDB/nr
log=nr.log

echo "nohup /data1/software/diamond-linux64/diamond blastx -q $pep -d $dbdir -o ./$out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &"
#10am-7pm

pep=../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa
out=uniprot.out
dbdir=/data1/Uniprot/uniprot_sprot
log=uniprot.log
#10min

#KOG
pep=../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa
out=kog.out
dbdir=/data1/KOGDatabase/clean.kog
log=kog.log

#real pep nr
pep=../1.Trinity/TransDecoder/Trinity_CD-HIT_0.9.fa.transdecoder.pep
out=nr_pep.out
dbdir=/data1/DiamondDB/nr
log=nr_pep.log

echo "nohup /data1/software/diamond-linux64/diamond blastp -q $pep -d $dbdir -o ./$out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &"


##################
##3.2 GO annotation & enrichment
##################

#eggnog注释: (GO KEGG COG/KOG)  (2.8W条序列--1hour)
emapper.py -i  $pep --output eggnog.out -m diamond --cpu 16 # 6.27 9am- 6.29 1am

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



##################
##3.3 KEGG anno & enrichment
##################

#backgroud
python /data1/script/GO_KEGG_Annotation/GetKEGGID_Annotation.py eggnog.out.emapper.annotations /data1/KEGG/ko00000_kegg.2018.5.30.txt > KEGG_Annotation.txt

cat KEGG_Annotation.txt|awk '
BEGIN{FS=OFS="\t"}
{split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4,$5}' > KEGG_Annotation_gene.txt

#sig genes
ln -s ../GO/sig.list ./

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$1;next}
{if($1 in a) print $0}' sig.list KEGG_Annotation_gene.txt > LvsA_KEGG_annotation.txt

#enrichment

python /data1/script/GO_KEGG_Annotation/Pyher_KEGGEnrichment_Analysis.py LvsA_KEGG_annotation.txt KEGG_Annotation_gene.txt LvsA


##################
##3.4 KAAS
##################
#KAAS
#1 chashuguzi
#2 freyasnow
#3 coconaut.freya
#4 iceyleemj

cat query.ko*.txt > ko.all.txt #14916 KAAS > eggnog 1475  

#只保留第二列不为空的
cat ko.all.txt|awk 'BEGIN{FS=OFS="\t"}$2!=""{print $0}' > temp

#将第二列输出为第7列。
cat temp|awk 'BEGIN{FS=OFS="\t"}{print $1,"","","","","",$2}'>ko.all4anno.txt

python /data1/script/GO_KEGG_Annotation/GetKEGGID_Annotation.py ko.all4anno.txt /data1/KEGG/ko00000_kegg.2018.5.30.txt > KEGG_Annotation.txt

#Enrichment
cat KEGG_Annotation.txt|awk '
BEGIN{FS=OFS="\t"}
{split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4,$5}' > KEGG_Annotation_gene.txt
ln -s ../GO_KEGG/sig.list ./

#保留_i1
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$1;next}
{split($1,x,"_");p=x[1]"_"x[2]"_"x[3]"_"x[4];if(p in a) print $0}' sig.list KEGG_Annotation.txt > LvsA_KEGG_annotation2.txt

python /data1/script/GO_KEGG_Annotation/Pyher_KEGGEnrichment_Analysis.py LvsA_KEGG_annotation2.txt KEGG_Annotation.txt LvsA


python /data1/script/GO_KEGG_Annotation/Pyher_KEGGEnrichment_Analysis.py LvsA_KEGG_annotation.txt KEGG_Annotation_gene.txt LvsA


##################
##3.4 Plot
##################
#自动识别当前目录下“.GO.Enrichment_out.txt” “.KEGG.Enrichment_out.txt”文件
Rscript /data1/script/GO_KEGG_Annotation/Plot_GO_KEGG.R "LvsA"



##################
##3.5 整合结果
##################
cd GO_KEGG
cat GO_annotation_gene.txt|cut -f 1,2| awk 'BEGIN { FS=OFS="\t" }
{
    curr = $1
    if (curr == prev) {
        rec = rec ";" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }' > Gene_GO.out.txt

cd KAAS
cat KEGG_Annotation_gene.txt|cut -f 1,2| awk 'BEGIN { FS=OFS="\t" }
{
    curr = $1
    if (curr == prev) {
        rec = rec ";" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }' > Gene_KEGG.out.txt



################################################
##4.Gene Expression
################################################

mkdir 4.Expression
mkdir Gene;mkdir Isoform

dir=../../3.Mapping

### Gene
cd Gene
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix RSEM --name_sample_by_basedir $dir/WR180739S_rsem_outdir/RSEM.genes.results $dir/WR180740S_rsem_outdir/RSEM.genes.results $dir/WR180741S_rsem_outdir/RSEM.genes.results $dir/WR180742S_rsem_outdir/RSEM.genes.results

### Isoform
cd Isoform
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix RSEM_isoform --name_sample_by_basedir $dir/WR180739S_rsem_outdir/RSEM.isoforms.results $dir/WR180740S_rsem_outdir/RSEM.isoforms.results $dir/WR180741S_rsem_outdir/RSEM.isoforms.results $dir/WR180742S_rsem_outdir/RSEM.isoforms.results 


################
##4.2 Diff Expression
################

###DESeq2 ### Diff Analysis use Gene 
### If you want see isoform,just use "RSEM_isoform.counts.matrix" 

cd 4.Expression/Gene
/data1/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.counts.matrix --method DESeq2 --samples_file sample_list.txt --output DESeq2_results --contrasts contrast.list

# L       WR180739S_rsem_outdir
# L       WR180740S_rsem_outdir
# A       WR180741S_rsem_outdir
# A       WR180742S_rsem_outdir

#contrast.list
#L	A

# L1	WR180739S 实验
# L2	WR180740S 实验
# A1	WR180741S 对照
# A2	WR180742S 对照

cd DESeq2_results/

#clustering  (差异基因 p<1e-3,fold>2^2)
/data1/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../RSEM.TMM.EXPR.matrix --samples ../sample_list.txt -P 1e-3 -C 2 --output cluster_results 

/data1/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../RSEM.TMM.EXPR.matrix --samples ../sample_list.txt -P 0.8 -C 0 --output cluster_results2 

#PCA
mkdir PCA;cd PCA
### PCA
/data1/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix ../RSEM.counts.matrix   -s ../sample_list.txt  --sample_cor_matrix  --prin_comp 3

################
##4.3 Deseq2
################

mkdir 5.Diff;cd 5.Diff

#RSEM.counts.matrix 修改首列，保证为样本名
#colData.csv 与matrix 列名样本顺序一致（最好为顺序）
cat RSEM.counts.matrix|awk 'BEGIN{FS="\t";OFS=","}{$1=$1;print $0}'>temp
#去掉小数点后两位
cat temp|sed "s/\...//g" >CountMatrix4DESeq.csv

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R CountMatrix4DESeq.csv colData.csv A ./


################################################
##5.stat & plot
################################################
cd 6.stats

### N50 
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl 1.Trinity/Trinity.fasta > 6.stats/Trinity_Stats.txt

/data1/software/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl 1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa > 6.stats/Trinity_CD-HIT_0.9_Stats.txt



#mapping stat
samtools flagstat myaln.bam

#FPKM
cd 4.Expression/Gene
#任意样本的结果提取基因长度信息
sample=WR180740S
cut -f 1,3,4 ../../3.Mapping/${sample}_rsem_outdir/RSEM.genes.results > $sample.feature_lengths.txt
#TMM方法将raw counts matrix转换为FPKM matrix
/data1/software/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix RSEM.counts.matrix --lengths $sample.feature_lengths.txt
#排序
cat RSEM.counts.matrix.TMM_normalized.FPKM|sort -k2,2 -n -r > FPKM.stats.txt


### Get Seq Length
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/misc/fasta_seq_length.pl ../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa  > Trinity_CD-HIT_0.9.seqLengths


################################################
##6.指定差异基因
################################################
cat nr.out|grep -E "[sS]ensory neuron membrane"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 > SNM.tsv
cat nr.out|grep -E "[iI]onotropic receptor"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 > IR.tsv
cat nr.out|grep -E "[cC]hemosensory protein"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 >CP.tsv
cat nr.out|grep -E "[gG]ustatory receptor"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 > GR.tsv
cat nr.out|grep -E "[oO]dorant binding protein"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 > OBP.tsv
cat nr.out|grep -E "[oO]lfactory receptor"|grep -v -E "Acyrthosiphon pisum|Sitobion avenae|Myzus persicae"|cut -f 1,2,13 > OR.tsv
#199






ln -s ../5.Diff/3.DiffExprGene/LvsA_all_genes_exprData.txt ./

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2"\t"$3;next}
{if ($1 in a) print $0,a[$1]}' SNMP_gene.list LvsA_all_genes_exprData.txt > SNMP_genes.txt

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2"\t"$3;next}
FNR==1{print $0,"GeneID","Description"}
FNR>1{if ($1 in a) print $0,a[$1]}' SNMP_gene.list LvsA_sig_genes_exprData.txt > SNMP_sig_genes.txt

ln -s ../4.Expression/Gene/RSEM.counts.matrix.TMM_normalized.FPKM ./
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$8"\t"$9;next}
FNR==1{print $0,"GeneID","Description"}
FNR>1{if ($1 in a) print $0,a[$1]}' SNMP_sig_genes.txt RSEM.counts.matrix.TMM_normalized.FPKM > SNMP_sig_FPKM.txt

ln -s ../5.Diff/CountMatrix4DESeq.csv ./



#按绿色和绿色+黄色筛选
for gene in SNM IR CP GR OBP OR
do
    awk '
    BEGIN{FS=OFS="\t"}
    NR==FNR{a[$1]=$1}
    {if ($3 in a)print $0}' ${gene}_G.list $gene.tsv > ${gene}_G.tsv
done

cat CountMatrix4DESeq.csv|awk 'BEGIN{FS=",";OFS="\t"}{$1=$1;print $0}' > temp1

for gene in SNM IR CP GR OBP OR
do
    cat ${gene}_G.tsv|awk '
    BEGIN{FS=OFS="\t"}
    {split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4}' > ${gene}_G_gene.list   
    
    # awk '
    # BEGIN{FS=OFS="\t"}
    # NR==FNR{a[$1]=$2;next}
    # FNR>1{if ($1 in a) print a[$1],$0}' ${gene}_G_gene.list  temp1|cut -f 1,3,4,5,6 > temp2

    # cat temp2|sort -k1,1 -u> temp3
    # cat temp3|awk 'BEGIN{FS="\t";OFS=",";print ",WR180739S,WR180740S,WR180741S,WR180742S"}{$1=$1;print $0}' >temp4

    # mv temp4 ${gene}_counts.csv


    #norm count matrix

    awk '
    BEGIN{FS=OFS="\t"}
    NR==FNR{a[$1]=$2;next}
    FNR>1{if ($1 in a) print a[$1],$0}' ${gene}_G_gene.list  norm-count-matrix.txt|cut -f 1,3,4,5,6> temp2

    cat temp2|sort -k1,1 -u|awk 'BEGIN{print "id\tWR180739S\tWR180740S\tWR180741S\tWR180742S"}{print $0}' > R_input/${gene}_G_normCountMatrix.tsv

    # mkdir $gene
    # Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R ${gene}_counts.csv colData.csv A ./$gene

done

for gene in SNM IR CP GR OBP OR
do
    Rscript heatmap.R ${gene}_G_normCountMatrix.tsv ${gene}_G_
done


###
for gene in IR CP GR OR
do
    awk '
    BEGIN{FS=OFS="\t"}
    NR==FNR{a[$1]=$1}
    {if ($3 in a)print $0}' ${gene}_GY.list $gene.tsv > ${gene}_GY.tsv
done

cat CountMatrix4DESeq.csv|awk 'BEGIN{FS=",";OFS="\t"}{$1=$1;print $0}' > temp1

for gene in IR CP GR OR
do
    cat ${gene}_GY.tsv|awk '
    BEGIN{FS=OFS="\t"}
    {split($1,x,"_");print x[1]"_"x[2]"_"x[3]"_"x[4],$2,$3,$4}' > ${gene}_GY_gene.list   

    awk '
    BEGIN{FS=OFS="\t"}
    NR==FNR{a[$1]=$2;next}
    FNR>1{if ($1 in a) print a[$1],$0}' ${gene}_GY_gene.list  norm-count-matrix.txt|cut -f 1,3,4,5,6> temp2

    cat temp2|sort -k1,1 -u|awk 'BEGIN{print "id\tWR180739S\tWR180740S\tWR180741S\tWR180742S"}{print $0}' > R_input/${gene}_GY_normCountMatrix.tsv
done

for gene in IR CP GR OR
do
    Rscript heatmap.R ${gene}_GY_normCountMatrix.tsv ${gene}_GY_
done

###2.Annotation/GO_KEGG 上下调分两列显示

###3.fasta annotation
#nr.out 13
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr.out  Trinity.fasta >Trinity_anno.fasta

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr_pep.out  Trinity_CD-HIT_0.9.fa.transdecoder.pep >Trinity_anno.pep.fasta










####gene name heatmap

cat gene_name.list |sort -k1,1 -u > gene_name_du.list
#去掉中括号
cat gene_name_du.list | sed "s/\[.*\]//"|sed "s/PREDICTED: //"|sed "s/LOW QUALITY PROTEIN: //" > temp
#保留
cat temp|awk '
BEGIN{FS=" ";OFS="\t"}{$1=$1}1'| awk '
BEGIN{FS=OFS="\t"}
{printf $1"\t";for(i=2;i<=NF;i++)printf$i" ";print ""}' >gene_name.tab

#大写改成小写
#删除refName字样

#conversion
for i in CP_G CP_GY GR_G GR_GY IR_G IR_GY OBP_G OR_G OR_GY SNM_G
do
  awk '
  BEGIN{FS=OFS="\t"}
  NR==FNR{a[$1]=$2;next}
  FNR>1{print a[$1],$2,$3,$4,$5}' gene_name.tab ${i}_normCountMatrix.tsv|sort -k1 > temp

  cat temp|sed "1 i Name\tWR180739S\tWR180740S\tWR180741S\tWR180742S" >nameTsv/${i}_name.tsv
done


for gene in SNM IR CP GR OBP OR
do
    Rscript heatmap.R ${gene}_G_name.tsv ${gene}_G_
done

for gene in IR CP GR OR
do
    Rscript heatmap.R ${gene}_GY_name.tsv ${gene}_GY_
done

