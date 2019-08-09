##############1.combine all hq&lq,and get fasta file##############
cat hq.fq >> lq.fq
mv lq.fq genome.fq
seqret -sequence genome.fq -outseq genome.fa

##############2.blastx: nuc vs prot##############
#nr,trembl,pfam,swissprot,go,kegg
#~/database/nr/nr cu06 12.21:16:00
#~/database/trembl/uniprot_trembl.fasta cu07 12.21:16:00
#~/database/pfam/Pfam-A-duplicate.fasta  cu05 12.21:16:00
#~/database/swissport/uniprot_sprot.fasta
#ls -sh -alt
cd analysis/pro4_isoseq/blast/
db=trembl
dbdir=~/database/trembl/uniprot_trembl.fasta
fasta=../genome.fa
log=${db}.log
out=${db}.out
echo $dbdir $fasta $log $out
nohup blastx -db $dbdir -query $fasta -num_threads 10 -out $out -outfmt "6" -evalue 1e-5 -max_target_seqs 1 > $log 2>&1 &

##############3.transdecoder##############
TransDecoder.LongOrfs -t ../genome.fa #30min 有进度百分比
TransDecoder.Predict -t ../genome.fa  #30-50min

cat genome.fa.transdecoder.cds|awk '
BEGIN{FS=OFS=" "}
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
!/>/{print}' > novel.cds.fa

cat genome.fa.transdecoder.pep|awk .. > novel.pep.fa

##############4.lncRNA##############
#get noncoding fasta :
#得到cds文件中所有isoform编号list
grep ">" genome.fa.transdecoder.cds |awk '{split($1,x,"::");print x[2]}' |sort -k1,1 -u > cds_isoform.list #108056 cds -> 53328 isoform
#在genome.fa去除cds.fa中的isoform
cat genome.fa |seqkit grep -f cds_isoform.list  -v >nc.fa #75415isoform - 53328 =22087

#与NR库比对
db=nr
dbdir=~/database/nr/nr
fasta=nc.fa
log=${db}.log
out=${db}.out
echo $dbdir $fasta $log $out
nohup blastx -db $dbdir -query $fasta -num_threads 10 -out $out -outfmt "6" -evalue 1e-5 -max_target_seqs 1 > $log 2>&1 &

#PELK (py2.7)
PLEKDIR=/home/liyb/biosoft/PLEK.1.2
python $PLEKDIR/PLEK.py -fasta nc.fa -out PLEK.out -thread 8 #Coding: 14363/22087=65.0292026984%, Non-coding: 7724/22087=34.9707973016%
cat PLEK.out |grep Non-coding >PLEK_nc.list  #7724
#cat PLEK.out |grep Coding|wc -l 14363
#cat PLEK.out |wc -l 22087


#CNCI(py2.7) 慢 5h
CNCIDIR=/home/liyb/biosoft/CNCI-master 
#-m 脊椎动物或者植物 参数 ve pl
python $CNCIDIR/CNCI.py -f nc.fa -o CNCI_out -p 10 
cat CNCI_out/CNCI.index |awk 'BEGIN{OFS="\t"}$2=="noncoding"{print $1,$2}' > CNCI_nc.list
#21173

#CPAT(py2.7)
make_hexamer_tab.py -c ../transdecoder/novel.cds.fa -n nc.fa > cpat.tsv
make_logitModel.py  -x cpat.tsv -c ../../transdecoder/novel.cds.fa -n ../nc.fa -o cpat
#?? 输入-g 是nc还是cds
cpat.py -g ../nc.fa  -d cpat.logit.RData -x cpat.tsv -o cpat_out
#score 小于0.36 认为non coding
cat cpat_out|awk 'BEGIN{OFS="\t"}NR>1&&$6<=0.36{print $1,$6}'> cpat_nc.list #22087->21291

#######Dr. Guo 采用给的模式生物做model（这里选择了mouse）
#物种  dat
species=Mouse
cpat.py -g ../../nc.fa -d ~/biosoft/CPAT-1.2.4/dat/${species}_logitModel.RData -x ~/biosoft/CPAT-1.2.4/dat/${species}_Hexamer.tsv -o CPAT_predicted



#合并处理

ln -s ../PLEK/PLEK_nc.list ./
ln -s ../CNCI_out/CNCI_nc.list ./
ln -s ../cpat/cpat_nc.list ./
#提取gene list
cat cpat_nc.list|awk '{split($1,x,"_");print tolower(x[1]) "_" x[2] "_" tolower(x[3])}' >cpat.list
cat CNCI_nc.list |awk 'BEGIN{FS=" "}{print $1}'> CNCI.list
cat PLEK_nc.list |awk '{split($2,x," ");print x[1]}'>temp
sed -i 's/>/ /g' temp
mv temp PLEK.list
#取交集
comm -12 <( sort CNCI.list ) <( sort PLEK.list ) >CNCI_PLEK.list #7572
comm -12 <( sort CNCI_PLEK.list ) <( sort cpat.list ) >CNCI_PLEK_cpat.list #6835
#提取lncRNA序列
cat ../../genome.fa |seqkit grep -f CNCI_PLEK_cpat.list > lncRNA.fasta


##############5.SSR##############
#MISA 输出到fasta所在目录：nc.fa 
MISADIR=/home/liyb/biosoft/MISA
cp $MISADIR/misa.ini ./
ln -s ../nc.fa ./
perl $MISADIR/misa.pl nc.fa







