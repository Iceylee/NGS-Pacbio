##############1.combine all hq&lq,and get fasta file##############
touch cDNA.fq
cat *.fq > cDNA.fq
seqret -sequence cDNA.fq -outseq cDNA.fa

cat cDNA.fa|grep ">"|wc -l #308298
cat cDNA.fq| grep -E "^@"|wc -l #质量值也有这个符号，因此计数不准确

cat *hq* > cDNA.hq.fq
seqret -sequence cDNA.hq.fq -outseq cDNA.hq.fa
seqret -sequence cDNA.lq.fq -outseq cDNA.lq.fa

#####cogent


##############2.diamond###########
cd /home/liyb/analysis/
db=trembl
fasta=../cDNA.unique.fa

dbdir=~/analysis/database/diamond_index/${db}.dmnd
log=$db.log
out=$db.out
echo $dbdir $fasta $log $out

nohup diamond blastx -d ${dbdir} -q $fasta -o $db.out -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive &

#--sensitive 

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

#CNCI(py2.7) 慢 5h
source activate Icey
CNCIDIR=/home/liyb/biosoft/CNCI-master 
#-m 脊椎动物或者植物 参数 ve pl
nohup python $CNCIDIR/CNCI.py -f nc.fa -o CNCI_out -p 15 &
cat CNCI_out/CNCI.index |awk 'BEGIN{OFS="\t"}$2=="noncoding"{print $1,$2}' > CNCI_nc.list

#CPAT(py2.7)
source activate Icey
species=Fly
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
mv na.fa.statistics SSR.statistics

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

#图片
cp -r /home/software/prog/smrtlink_data/jobs/000/000317/html/images ~/analysis/pro10_maie_isoseq/images/90
cp -r /home/software/prog/smrtlink_data/jobs/000/000318/html/images ~/analysis/pro10_maie_isoseq/images/91

cd 90
mv images/pbreports.tasks.*/* /
rm -f *thumb*
for i in `ls `
do
	mv $i under4k_$i
done

cd 91
mv images/pbreports.tasks.*/* /
rm -f *thumb*
for i in `ls `
do
	mv $i above4k_$i
done

#stats
cd Diamond/
wc -l *.anno.out|grep -o -E "[0-9]+"|xargs -I {} echo {}/50334|bc -l






