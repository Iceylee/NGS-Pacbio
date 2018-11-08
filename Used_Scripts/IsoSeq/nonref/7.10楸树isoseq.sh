

##############1.combine all hq&lq,and get fasta file##############
cd raw/
touch cDNA.fq
cat *.fastq > cDNA.fq
seqret -sequence cDNA.fq -outseq cDNA.fa

cat cDNA.fa|grep ">"|wc -l ##cDNA.fa 319917

##############2.cogent##############
17      5289
44      7608

wc -l cmd  #57564/3/6 = 序列数 （若小数加1凑整）*3
split -l 9600 cmd -d

cd analysis/pro16_tree/cogent
source activate anaCogent
#分节点运行x00 到 x05
#时刻查看top 如果某个任务运行时间过长（且在>4M 记录中） kill 记录bin的编号
ps -awx|grep 12332  
#kill:5289 10219 7608



#kill的bin 再次分bin
num=1 #后续不改，延续文件夹编号

max=40000
parts=6

for folder in 5289
do
	echo "cd $folder"
	echo "perl ~/biosoft/fasta-splitter.pl --n-parts $parts isoseq_flnc.fasta" 
	echo "cd .."
	for i in $(seq 1 $parts) 
	do
	Fold=`expr $num + $max` 
	echo "mkdir $Fold"
	echo "mv $folder/isoseq_flnc.part-${i}.fasta $Fold/isoseq_flnc.fasta"
	num=`expr $num + 1`
	done
done

parts=12
for folder in 7608
do
	echo "cd $folder"
	echo "perl ~/biosoft/fasta-splitter.pl --n-parts $parts isoseq_flnc.fasta" 
	echo "cd .."
	for i in $(seq 1 $parts) 
	do
	Fold=`expr $num + $max` 
	echo "mkdir $Fold"
	echo "mv $folder/isoseq_flnc.part-${i}.fasta $Fold/isoseq_flnc.fasta"
	num=`expr $num + 1`
	done
done

#5289 40001-4006
#7608 40007 -40018


bash ~/analysis/scripts/Cogent_check1.sh 

#去除原始bin
mv preCluster_out/5289/ temp/
mv preCluster_out/7608/ temp/

#cmd2失败的重跑
cat cmd2.failed.list |awk -F '[" "//]' '{printf "reconstruct_contig.py cDNA/"$2" -p "$2"\n"}' >cmd2.last #316
split -l 4 ../cmd2.last -d

#再失败的 添加参数k 50 
bash ~/analysis/scripts/Cogent_check2.sh  > cmd2.failed.list
cat cmd2.failed.list |grep -E -o "[0-9]+_[0-9]+" > cmd2.failed.bin

for line in `cat cmd2.failed.bin`
do
 echo "nohup reconstruct_contig.py cDNA/$line -p $line -k 50 &"
done

#查看hello.log 调整k值（慢慢变大从40） 不一定50

#AssertionError #认为polyA尾没有去掉 旧的cogent版本
#去了发现并没有polyA尾巴
https://github.com/bowhan/trim_isoseq_polyA.git


# mkdir build && \
# cd build && \
# cmake ../ -DBOOST_ROOT=/home/software/soft/boost_1_63_0 -DCMAKE_BUILD_TYPE=Release && \
# make

# /home/liyb/biosoft/trim_isoseq_polyA/bin/trim_isoseq_polyA -i cDNA/10266_0/in.fa -t 8 -G > cDNA/10266_0/input.atrim.fa 2> cDNA/10266_0/input.atrim.log


#16failed

##############3.diamond###########
#nr,trembl,pfam
#sp,cog

cd /home/liyb/analysis/pro16_tree/Diamond
db=sp
fasta=../cDNA.unique.fa

dbdir=~/analysis/database/diamond_index/${db}.dmnd
log=$db.log
out=$db.out
echo $dbdir $fasta $log $out

nohup diamond blastx -d ${dbdir} -q $fasta -o $db.out -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &

#--sensitive 

for db in nr trembl pfam sp cog
do
	cat $db.out|sort -k1,1 -u > $db.anno.out

done

mkdir result
mv *.anno.out result/

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

#######anno fasta
ln -s ../Diamond/nr.out ./
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr.out  ../cDNA.unique.fa > cDNA_anno.fasta

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr.out  novel.pep.fa >pep_anno.fasta

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' nr.out  novel.cds.fa >cds_anno.fasta

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

#合并
mv *nc.list result/
cd result
#提取gene list
#不同文库的情况

cat PLEK_nc.list |awk '{split($3,x," ");print x[1]}'|sed 's/>//g' >PLEK.list
cat CNCI_nc.list |awk 'BEGIN{FS=" "}{print $1}'> CNCI.list

#取交集
comm -12 <( sort CNCI.list ) <( sort PLEK.list ) >CNCI_PLEK.list 

##去掉比对上pfam库的序列
ln -s ../../Diamond/pfam.out
cat pfam.out|cut -f 1|sort -k1,1 -u > pfam.list
#仅在CNCI_PLEK_CPAT.list 中出现的序列
comm -23 <( sort CNCI_PLEK.list ) <( sort pfam.list ) >CNCI_PLEK_pfam.list 

#提取lncRNA序列
cat ../../cDNA.unique.fa |seqkit grep -f CNCI_PLEK_pfam.list  > lncRNA.original.fasta

#ORF预测 有编码潜能的序列
mkdir angel && cd angel
source activate Icey
dumb_predict.py ../lncRNA.original.fasta angel.dumb --cpus 15
#--min_aa_length 300 

cat angel.dumb.final.cds|grep ">" |cut -d " " -f 1|sed 's/>//g'|sed 's/|m.*//g' > angel.list
#去掉angel.list中的序列(23条)
cat ../lncRNA.original.fasta |seqkit grep -f angel.list -v > lncRNA.fasta



##############5.SSR##############
#MISA 输出到fasta所在目录：nc.fa 
ln -s ../lncRNA/nc.fa
MISADIR=/home/liyb/biosoft/MISA
cp $MISADIR/misa.ini ./
nohup perl $MISADIR/misa.pl nc.fa &
mv nc.fa.misa SSR.misa
mv nc.fa.statistics SSR.statistics

#########6.pfam18##############
cd pfam18
ln -s ../Diamond/result/pfam.anno.out
ln -s ~/analysis/database/pfam18
cat pfam.anno.out|cut -f 1,2,11 > temp

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$2]=$4;next;}
{if (a[$2]!="") print $0,a[$2]}
' pfam18 temp > temp2

sed '1 i Sequence_ID\tHMM_ID\tE_value\tDescription' temp2 > pfam18.annotation.txt






#########results##############
##
cat cDNA.unique.fa|grep ">" > temp.list
sed -i 's/>//g' temp.list #在原文件删除 匹配大于号删除
cat temp.list|awk 'BEGIN{FS=" "}{print $1}' > genome.list 

#图片
#分析项目目录下
cp -r /opt/pacbio_soft/smrtlink5.0/smrtlink/userdata/jobs_root/000/000409/html/images ./images/



mv images/pbreports.tasks.*/* ./
rm -f *thumb*
for i in `ls *.png`
do
	mv $i under4k_$i
done




# stats
cd Diamond/
wc -l *.anno.out
wc -l *.anno.out|grep -o -E "[0-9]+"|xargs -I {} echo {}/54280|bc -l

wc -l pfam18.annotation.txt
 



