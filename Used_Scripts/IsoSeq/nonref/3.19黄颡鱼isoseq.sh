##############1.combine all hq&lq,and get fasta file##############
cd raw/
touch cDNA.fq
cat *.fastq > cDNA.fq
seqret -sequence cDNA.fq -outseq cDNA.fa

cat cDNA.fa|grep ">"|wc -l #452020

##############2.cogent##############

##############3.diamond###########
#nr,trembl,pfam
#sp,cog

cd ~/analysis/pro11_yellowFish_isoseq/Diamond
db=cog
fasta=../cDNA.unique.fa

dbdir=~/analysis/database/diamond_index/${db}.dmnd
log=$db.log
out=$db.out
echo $dbdir $fasta $log $out

nohup diamond blastx -d ${dbdir} -q $fasta -o $db.out -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive &

#--sensitive 

#diamond blastx -q ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa -d /data/KOGDatabase/clean.kog  -o KOG_annotation.bk.out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle


##############3.transdecoder##############
mkdir TransDecoder
cd TransDecoder

TransDecoder.LongOrfs -t ../cDNA.unique.fa #10min
TransDecoder.Predict -t ../cDNA.unique.fa  #30-50min




# #sed 元字符需要反义 grep 非元字符需要翻译
# cat cDNA.unique.fa.transdecoder.cds|sed 's/|[0-9]\+_[0-9]\+|path[0-9]\+:[0-9]\+-[0-9]\+(.)//g' > del.cds
# cat cDNA.unique.fa.transdecoder.pep|sed 's/|[0-9]\+_[0-9]\+|path[0-9]\+:[0-9]\+-[0-9]\+(.)//g' > del.pep



echo 'BEGIN{FS=OFS=" "}
/>/{split($7,x,":");
s1=">" x[1];
s2=$4;
s3=$5;
gsub("[\\(\\)]","",$6);
s4="strand:" $6;
split(x[2],y,"(");
s5="pos:" y[1]
print s1,s2,s3,s4,s5
}
!/>/{print}'  > name.awk

cat cDNA.unique.fa.transdecoder.cds|awk -f name.awk > novel.cds.fa
cat cDNA.unique.fa.transdecoder.pep|awk -f name.awk > novel.pep.fa



##############4.lncRNA##############
ln -s ../TransDecoder/novel.cds.fa
# get noncoding fasta :
# 得到cds文件中所有isoform编号list
grep ">" novel.cds.fa|cut -d " " -f 1|sed 's/>//g'|sort -k1,1 -u > cds_isoform.list #138218
# 在cDNA.fa去除cds.fa中的isoform
cat ../cDNA.unique.fa |seqkit grep -f cds_isoform.list  -v > nc.fa 

#CPC很慢 

#pfamscan


#与NR库比对？如何确定lncRNA与否

#PELK (py2.7)
source activate Icey
PLEKDIR=/home/liyb/biosoft/PLEK.1.2
nohup python $PLEKDIR/PLEK.py -fasta nc.fa -out PLEK.out -thread 15 1>PLEK.log &
cat PLEK.out |grep Non-coding >PLEK_nc.list  
#exit 但log显示complete

#CNCI(py2.7) 慢 5h
source activate Icey
CNCIDIR=/home/liyb/biosoft/CNCI-master 
#-m 脊椎动物或者植物 参数 ve pl
nohup python $CNCIDIR/CNCI.py -f nc.fa -o CNCI_out -p 15 &
cat CNCI_out/CNCI.index |awk 'BEGIN{OFS="\t"}$2=="noncoding"{print $1,$2}' > CNCI_nc.list

#CPAT(py2.7)
source activate Icey
species=Zebrafish #Zebrafish Human Mouse Fly
#0.38 0.364 0.44 0.39

nohup cpat.py -g nc.fa -d ~/biosoft/CPAT-1.2.4/dat/${species}_logitModel.RData -x ~/biosoft/CPAT-1.2.4/dat/${species}_Hexamer.tsv -o CPAT.out &
cat CPAT.out|awk 'BEGIN{OFS="\t"}NR>1&&$6<=0.38{print $1,$6}'> CPAT_nc.list 

#合并
mv *nc.list result/
#提取gene list
#不同文库的情况

cat CPAT_nc.list|awk '
{split($1,x,"|");
	split(x[2],y,"_");
	sample = tolower(y[1]) "_" y[2] "_" tolower(y[3])
	print(x[1]"|"sample"|"tolower(x[3]))
}' >CPAT.list
#cat CPAT_nc.list|awk '{print tolower($1)}' > CPAT.list
cat PLEK_nc.list |awk '{split($3,x," ");print x[1]}'|sed 's/>//g' >PLEK.list
cat CNCI_nc.list |awk 'BEGIN{FS=" "}{print $1}'> CNCI.list

#取交集
comm -12 <( sort CNCI.list ) <( sort PLEK.list ) >CNCI_PLEK.list 
comm -12 <( sort CNCI_PLEK.list ) <( sort CPAT.list ) >CNCI_PLEK_CPAT.list 


##去掉比对上pfam库的序列
ln -s /home/liyb/analysis/pro10_maie_isoseq/Diamond/result/pfam.anno.out
cat pfam.anno.out|cut -f 1 > pfam.list
#仅在CNCI_PLEK_CPAT.list 中出现的序列
comm -23 <( sort CNCI_PLEK_CPAT.list ) <( sort pfam.list ) >CNCI_PLEK_CPAT_pfam.list 

#提取lncRNA序列
cat ../../cDNA.unique.fa |seqkit grep -f CNCI_PLEK_CPAT_pfam.list > lncRNA.original.fasta

#ORF预测 有编码潜能的序列
dumb_predict.py lncRNA.original.fasta angel.dumb --cpus 20 #--min_aa_length 300 
cat angel.dumb.final.cds|grep ">" |cut -d " " -f 1|sed 's/>//g'|sed 's/|m.*//g' > angel.list
#去掉angel.list中的序列 
cat lncRNA.original.fasta |seqkit grep -f angel.list -v > lncRNA.fasta


##############5.SSR##############
#MISA 输出到fasta所在目录：nc.fa 
ln -s ../lncRNA/nc.fa
MISADIR=/home/liyb/biosoft/MISA
cp $MISADIR/misa.ini ./
nohup perl $MISADIR/misa.pl nc.fa &
mv nc.fa.misa SSR.misa
mv nc.fa.statistics SSR.statistics

#########6.pfam18
cd pfam18
ln -s ../Diamond/result/pfam.anno.out
ln -s ~/analysis/database/pfam18
cat pfam.anno.out|cut -f 1,11,13|awk '
BEGIN{FS=OFS="\t"}
{
split($3,x," ")
if (x[1] != "") {print $1,x[1],$2}
}' > temp

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$4;next;}
{if (a[$2]!="") print $0,a[$2]}
' pfam18 temp > temp2

sed '1 i Sequence_ID\tHMM_ID\tE_value\tDescription' temp2 > pfam18.annotation.txt

##
cat cDNA.unique.fa|grep ">" > temp.list
sed -i 's/>//g' temp.list #在原文件删除 匹配大于号删除
cat temp.list|awk 'BEGIN{FS=" "}{print $1}' > genome.list 

#图片
#分析项目目录下
cp -r /home/software/prog/smrtlink_data/jobs/000/000334/html/images ./images/178
cp -r /home/software/prog/smrtlink_data/jobs/000/000336/html/images ./images/179


cp 178/pbreports.tasks.*/* ./
rm -f *thumb*
for i in `ls *.png`
do
	mv $i under4k_$i
done


cp 179/pbreports.tasks.*/* ./
rm -f *thumb*
for i in `ls *.png|grep -v "under4k"`
do
	mv $i above4k_$i
done


# stats
cd Diamond/
wc -l *.anno.out
wc -l *.anno.out|grep -o -E "[0-9]+"|xargs -I {} echo {}/97609|bc -l
wc -l pfam.anno.out #-1

 
## 有参
cd gmap
nohup gmap -D /home/liyb/miniconda3/share/Danio_rerio -d Danio_rerio -f samse -n 0 -t 20 ../raw/cDNA.fa > gmap.sam 2> gmap.sam.log &

cd gmap2
nohup gmap -D /home/liyb/miniconda3/share/Danio_rerio -d Danio_rerio -f samse -n 0 -t 20 ../cDNA.unique.fa > gmap.sam 2> gmap.sam.log &
#分开重跑
nohup gmap -D /home/liyb/miniconda3/share/Danio_rerio -d Danio_rerio -f samse -n 0 -t 20 ../2.fa > gmap.2.sam 2> gmap.2.sam.log &


# 斑点叉尾鮰 下载genome fasta（不是转录本）
mv GCF_001660625.1_IpCoco_1.2_genomic.fna Ictalurus_punctatus.fa
gmap_build -d Ictalurus_punctatus raw/Ictalurus_punctatus.fa

cd gmap #cu07
nohup gmap -D /home/liyb/miniconda3/share/Ictalurus_punctatus -d Ictalurus_punctatus -f samse -n 0 -t 20 ../raw/cDNA.fa > gmap.Ictalurus_punctatus.sam 2> gmap.Ictalurus_punctatus.sam.log &

cd gmap2 #cu08
nohup gmap -D /home/liyb/miniconda3/share/Ictalurus_punctatus -d Ictalurus_punctatus -f samse -n 0 -t 20 ../cDNA.unique.fa > gmap.Ictalurus_punctatus.sam 2> gmap.Ictalurus_punctatus.sam.log &
#重跑2
nohup gmap -D /home/liyb/miniconda3/share/Ictalurus_punctatus -d Ictalurus_punctatus -f samse -n 0 -t 20 ../2.fa > gmap.Ictalurus_punctatus.2.sam 2> gmap.Ictalurus_punctatus.2.sam.log &

###
cat gmap.all.sam |sort -k1,1 -u|wc -l
#93563
cat gmap.all.sam |sort -k1,1 -u|cut -f 2|grep -v 4|wc -l
#50631

cat gmap.Ictalurus_punctatus.all.sam |sort -k1,1 -u|wc -l
#97612
cat gmap.Ictalurus_punctatus.all.sam |sort -k1,1 -u|cut -f 2|grep -v 4|wc -l
#90152


##STAR
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir /home/liyb/ref/index/STAR/Ictalurus_punctatus --genomeFastaFiles ../raw/Ictalurus_punctatus.fa  --sjdbGTFfile ../raw/Ictalurus_punctatus.gff --sjdbOverhang 100 --sjdbGTFtagExonParentTranscript Parent
# --runThreadN ：设置线程数
# --genomeDir ：index输出的路径
# --genomeFastaFiles ：参考基因组序列
# --sjdbGTFfile ：参考基因组注释文件
# --sjdbOverhang ：这个是reads长度的最大值减1，默认是100 #pacbio也用100就可以了
# --sjdbGTFtagExonParentTranscript :gff 而不是gtf的时候需要此参数


STARlong --runMode alignReads \
--runThreadN 20 --genomeDir /home/liyb/ref/index/STAR/Ictalurus_punctatus \
--readFilesIn ../cDNA.unique.fa \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000
# --readFilesIn ：paired reads文件
# --outSAMtype ：表示输出默认排序的bam文件，类似于samtools sort（还有–outSAMtype BAM Unsorted和–outSAMtype BAM Unsorted SortedByCoordinate）
# --outFileNamePrefix ：输出文件路径即前缀


samtools view SRR3589959Aligned.sortedByCoord.out.bam |less -S

 
STARlong --runMode alignReads \
--runThreadN 20 --genomeDir /home/liyb/ref/index/STAR/Ictalurus_punctatus \
--readFilesIn ../raw/cDNA.fa \
--outSAMtype BAM SortedByCoordinate \
--outSAMattributes NH HI NM MD \
--readNameSeparator space \
--outFilterMultimapScoreRange 1 \
--outFilterMismatchNmax 2000 \
--scoreGapNoncan -20 \
--scoreGapGCAG -4 \
--scoreGapATAC -8 \
--scoreDelOpen -1 \
--scoreDelBase -1 \
--scoreInsOpen -1 \
--scoreInsBase -1 \
--alignEndsType Local \
--seedSearchStartLmax 50 \
--seedPerReadNmax 100000 \
--seedPerWindowNmax 1000 \
--alignTranscriptsPerReadNmax 100000 \
--alignTranscriptsPerWindowNmax 10000




