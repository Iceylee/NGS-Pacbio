#!/bin/bash

Num=1
for i in $(seq 1 4)
do
	if [ $Num -le 2 ];then
		Num=`expr $Num + 1`
		# cd /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}
		# run_mash.py -k 30 --cpus=25 /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta
		# process_kmer_to_graph.py /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta.s1000k30.dist /home/liyb/analysis/pro10_maie_isoseq/cogent/cDNA ${i}
		#echo ${i}
	else
		wait
		# cd /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}
		# run_mash.py -k 30 --cpus=25 /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta
		# process_kmer_to_graph.py /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta /home/liyb/analysis/pro10_maie_isoseq/cogent/preCluster_out/${i}/isoseq_flnc.fasta.s1000k30.dist /home/liyb/analysis/pro10_maie_isoseq/cogent/cDNA ${i}
		Num=2
	fi
		
done

