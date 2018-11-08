for i in 01 02 03 04 05 06 07 08 09 $(seq 10 40)
do
Fold=`expr $i + 20024`
echo "mkdir $Fold"
echo "mv 11917/isoseq_flnc.part-${i}.fasta $Fold/isoseq_flnc.fasta"
done

for i in $(seq 1 3)
do
Fold=`expr $i + 20064`
echo "mkdir $Fold"
echo "mv 7913/isoseq_flnc.part-${i}.fasta $Fold/isoseq_flnc.fasta"
done

#!/bin/bash
for i in $(seq 20065 20067)
do

	echo "cd /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}"
	echo "run_mash.py -k 30 --cpus=25 /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta"
	echo "process_kmer_to_graph.py /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta.s1000k30.dist /home/liyb/analysis/pro10_maie_isoseq/cogent/cDNA ${i}"
done


for i in $(seq 20005 20024) 11917 7913 8048
do
	echo "mv $i ~/temp/"
done

# Num=1
# for i in `ls cDNA`
# do	
# 	if Num/2 != 0 then
# 		reconstruct_contig.py cDNA/$i -p $i &
# 		Num=`expr $Num + 1`
# 	else
# 		wait
# 		reconstruct_contig.py cDNA/$i -p $i &
# done


# if [ ! -d "$folder"]; then


##cmd2
split -l 100 cmd2 -d
for i in $(seq 9011 9024)
do
	echo "nohup bash x${i} &"
done

