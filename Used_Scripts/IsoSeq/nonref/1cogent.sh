

################1.cluster#############
#run_preCluster.py --cpus=25   #3:51
run_preCluster.py --cpus=25 --dun_use_partial
generate_batch_cmd_for_Cogent_family_finding.py --cpus=25 --cmd_filename=cmd preCluster.cluster_info.csv preCluster_out cDNA
bash cmd 

################2.reconstruction#############
printf "Partition\tSize\tMembers\n" > final.partition.txt
ls preCluster_out/*/*partition.txt | xargs -n1 -i sed '1d; $d' {} | cat >> final.partition.txt
generate_batch_cmd_for_Cogent_reconstruction.py cDNA/ > cmd2
bash cmd2

# split -l 5466 cmd






################4.fake genome#############
export PATH=$PATH:~/biosoft/cDNA_Cupcake-master/sequence/
tail -n 1 final.partition.txt | tr ',' '\n'  > unassigned.list
get_seqs_from_list.py isoseq_flnc.fasta unassigned.list > unassigned.fasta

mkdir collected
cd collected
cat ../cDNA/*/cogent2.renamed.fasta ../unassigned.fasta > cogent.fake_genome.fasta


################5.gmap & collapse#############
# gmap
gmap_build -D . -d fake_genome cogent.fake_genome.fasta #1min 
gmap -D . -d fake_genome -f samse -n 0 -t 22 \
   ../isoseq_flnc.fasta > isoseq_flnc.fasta.sam 2> \
   isoseq_flnc.fasta.sam.log  #5min
#collapse
sort -k 3,3 -k 4,4n isoseq_flnc.fasta.sam > isoseq_flnc.fasta.sorted.sam
collapse_isoforms_by_sam.py --input ../isoseq_flnc.fasta -s isoseq_flnc.fasta.sorted.sam \
           -c 0.95 -i 0.85 --dun-merge-5-shorter \
           -o isoseq_flnc.fasta.no5merge
#没有cluster report 以下两步省略（非必须）
#下载cluster report并合并 
ln -s ../../raw/cluster_report.csv
#报错 
get_abundance_post_collapse.py isoseq_flnc.fasta.no5merge.collapsed cluster_report.csv
filter_away_subset.py isoseq_flnc.fasta.no5merge.collapsed

#isoseq_flnc.fasta.no5merge.collapsed.rep.fa 即输出fasta
cd ../../
ln -s hqCogent/collected/isoseq_flnc.fasta.no5merge.collapsed.rep.fa cDNA.unique.fa



####cDNA check1.sh
#!/bin/bash

counter=0
for d in preCluster_out/* ; do
    arr=(${d/\// })
    FILE=$d/${arr[1]}".partition.txt";
    if [ ! -f $FILE ]; then
        echo "$d failed. Please re-run."
        counter=$((counter+1))
    fi
done

if [ $counter -eq 0 ]; then
    echo "All bins completed! Proceed to reconstruction."
else
    echo "Some bins failed. Please re-run them."
fi


#check2.sh
#!/bin/bash

counter=0
DIRNAME='cDNA'

for d in $DIRNAME/* ; do
    FILE=$d/cogent2.fa
    if [ ! -f $FILE ]; then
        echo "$d failed. Please re-run."
        counter=$((counter+1))
    fi
done

if [ $counter -eq 0 ]; then
    echo "All jobs completed! Proceed to reconstruction."
else
    echo "Some jobs failed. Please re-run them."
fi





for i in $(seq 20005 20024) 11917 8048 
