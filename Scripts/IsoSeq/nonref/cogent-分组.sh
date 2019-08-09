for i in 01 02 03 04 05 06 07 08 09 $(seq 10 40)
do
Fold=`expr $i + 20024`
echo "mkdir $Fold"
echo "mv 11917/isoseq_flnc.part-${i}.fasta $Fold/isoseq_flnc.fasta"
done







#!/bin/bash
for i in $(seq 20065 20067)
do

	echo "cd /home/liyb/analysis/pro10_maie_isoseq/cogent/pnohup reCluster_out/${i}" &
	echo "run_mash.py -k 30 --cpus=25 /home/liyb/analysis/pro10_maie_isoseq/cogent/pnohup reCluster_out/${i}/ &isoseq_flnc.fasta"
	echo "process_kmer_to_graph.py /home/liyb/analysis/pro10_maie_isoseq/cogent/pnohup reCluster_out/${i}/isoseq_flnc.fasta /home/liyb/ &analysis/pro10_maie_isoseq/cogent/pnohup reCluster_out/${i}/ &isoseq_flnc.fasta.s1000k30.dist /home/liyb/analysis/pro10_maie_isoseq/cogent/cDNA ${i}"
done


for i in $(seq 20005 20024) 11917 7913 8048
do
	echo "mv $i ~/temp/"
done


####3.19
#pnohup recluster 32463 &


######1. 太大的fasta切分
du -sm * | sort -n > ../folder.size #(大于4M则分开)
ls |sort -n -k1,1|less #确定bin的最大编号

num=1 #后续不改，延续文件夹编号

max=40000
parts=6

for folder in 11805 19663 2133 4167 5894 15298 10219 19990 5289 7608 
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


###2.cmd生成 (只生成自动分的bin的cmd)
generate_batch_cmd_for_Cogent_family_finding.py --cpus=20 --cmd_filename=cmd preCluster.cluster_info.csv preCluster_out cDNA &
#分成6个
wc -l cmd  #77610/3/6 = 序列数 （若小数加1凑整）*3
split -l 12936 cmd -d

###3.自己分的bin 生成cmd
path=/home/liyb/analysis/pro16_tree/
for i in $(seq 40001 40018)
do
	echo "cd ${path}/cogent/preCluster_out/${i}" 
	echo "run_mash.py -k 30 --cpus=25 ${path}/cogent/preCluster_out/${i}/isoseq_flnc.fasta" 
	echo "process_kmer_to_graph.py /${path}/cogent/preCluster_out/${i}/isoseq_flnc.fasta ${path}/cogent/preCluster_out/${i}/isoseq_flnc.fasta.s1000k30.dist ${path}/cogent/cDNA ${i}"
done > cmd1_last

split -l 18 cmd1_last -d

for i in $(seq 0 4)
do
	echo "nohup bash x0${i} &"
done


#去掉split的那些bin的原始bin
for i in 29099 24665 21807 40001 40002 $(seq 40005 40016) $(seq 40020 40022)
do
	echo "mv ${i} ../temp/"
done



##cmd2

## 20 runs per cu
## 分成120组 19382/120=162
split -l 162 ../cmd2 -d

for i in $(seq 60 79)
do
	echo "nohup bash cmd2_last/x${i} &"
done

for i in $(seq 9010 9029 )
do
  echo "nohup bash cmd2_dir/x${i} &"
done

#将40000开头的提出来单独跑

nohup reconstruct_contig.py cDNA/40003_0 -p 40003_0 &
nohup reconstruct_contig.py cDNA/40003_2883 -p 40003_2883 &
nohup reconstruct_contig.py cDNA/40003_510 -p 40003_510 &
nohup reconstruct_contig.py cDNA/40003_972 -p 40003_972 &
nohup reconstruct_contig.py cDNA/40004_0 -p 40004_0 &
nohup reconstruct_contig.py cDNA/40004_3305 -p 40004_3305 &
nohup reconstruct_contig.py cDNA/400_0 -p 400_0 &
nohup reconstruct_contig.py cDNA/40017_0 -p 40017_0 &
nohup reconstruct_contig.py cDNA/40017_400 -p 40017_400 &
nohup reconstruct_contig.py cDNA/40018_0 -p 40018_0 &
nohup reconstruct_contig.py cDNA/40019_0 -p 40019_0 &
nohup reconstruct_contig.py cDNA/40019_2131 -p 40019_2131 &
nohup reconstruct_contig.py cDNA/4002_0 -p 4002_0 &
nohup reconstruct_contig.py cDNA/40023_0 -p 40023_0 &
nohup reconstruct_contig.py cDNA/40024_0 -p 40024_0 &
nohup reconstruct_contig.py cDNA/40025_0 -p 40025_0 &
nohup reconstruct_contig.py cDNA/40026_0 -p 40026_0 &
nohup reconstruct_contig.py cDNA/40027_0 -p 40027_0 &
nohup reconstruct_contig.py cDNA/40028_0 -p 40028_0 &


for i in cog go nr pfam trembl
do
	echo "mv ${i}.anno.out temp/"
	echo "mv ${i}.temp ${i}.anno.out"
done


#k 50

cat cmd2.failed.list |grep -E -o "[0-9]+_[0-9]+" > cmd2.failed.bin

for line in `cat cmd2.failed.bin`
do
 echo "nohup reconstruct_contig.py cDNA/$line -p $line -k 50 &"
done

#查看hello.log 是否需要挑高k 
#查看in.fa 是否有polyA

nohup reconstruct_contig.py cDNA/10641_5 -p 10641_5 -k 50 &
#max(arg)




nohup reconstruct_contig.py cDNA/11810_2 -p 11810_2 -k 50 &
#list index out of range
/home/liyb/biosoft/trim_isoseq_polyA/bin/trim_isoseq_polyA -i in.fa -t 8 -G > input.atrim.fa 2> input.atrim.log
mv in.fa temp.fa
mv input.atrim.fa in.fa
reconstruct_contig.py --nx_cycle_detection -k 50 . 




cDNA/13840_3
#max() arg is an empty sequence
cDNA/14485_0
#单个序列 k10 到50 都找不到
cDNA/18577_4
cDNA/19897_0
cDNA/21444_2
cDNA/21483_2
cDNA/2190_5
cDNA/22175_0
cDNA/22507_0
cDNA/6744_7
cDNA/12727_2
cDNA/13667_0 


#大序列
cDNA/15298_0
reconstruct_contig.py --nx_cycle_detection -k 70 . 

cDNA/40018_0
reconstruct_contig.py --nx_cycle_detection -k 60 . 

cDNA/5894_0
reconstruct_contig.py --nx_cycle_detection -k 60 . 





