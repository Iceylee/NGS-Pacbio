
########mummer##########

##1.将参考基因组和query对齐起点
nucmer --maxmatch -c 90 -l 40 ../other/KT2440_genome.fna  ../assembly.fasta 
show-coords -rclT out.delta > ref_qry.coords

#由coords，发现反向,调整
revseq A316.fasta A316.rev.fasta 
nucmer --maxmatch -c 90 -l 40 ../other/KT2440_genome.fna  ../assembly.fasta 
show-coords -rclT out.delta > ref_qry.coords

#按coords文件，拼接reference genome
#5.7M-3240130，5.7M +  1，5.7M-3240139
cat A316.rev.fasta |seqret --filter -sbegin 1 -send 2475684 > 1.fasta
cat A316.rev.fasta |seqret --filter -sbegin 2475685 -send 5715815 > 2.fasta
cat 2.fasta 1.fasta |union -filter > A316_paste.fasta


##2.plot方法-wechat
nucmer --maxmatch -c 90 -l 40 ../other/A316_paste.fasta  ../assembly.fasta 
#共线性作图
mummerplot out.delta -R ../other/A316_paste.fasta -Q ../assembly.fasta --png  --filter --layout


#3.plot方法-Dr. Guo + SNP
#输出文件 out.delta
nucmer --maxmatch -c 90 -l 40 ../other/A316_paste.fasta ../assembly.fasta   
#SNP 来自out.delta
delta-filter -r -q out.delta > out-filter.delta 

#plot 
mummer -mum -b -c -l 30 ../other/A316_paste.fasta ../assembly.fasta > out.mummer 
mummerplot -postscript out.mummer #out.ps
ps2pdf out.ps out.pdf

#snp&Indel
# -C option, only SNPs from uniquely aligned regions will be reported.
#show-snps out-filter.delta > SNP_Indel.txt
show-snps -Clr -x 1 -T out-filter.delta  >SNP_Indel.txt
#统计indel个数（DIST为1且P2不变。且总长<50)

########mauve##########
#文件路径为英文

####gblocks
Gblocks allSingleCopyOrthologsAlign.fasta -t=p

#mysql
mysql server start

