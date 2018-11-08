
################1.cluster#############
cd cogent
source activate anaCogent
#序列数>20,000
ln -s ../raw/cDNA.fa isoseq_flnc.fasta
nohup run_preCluster.py --cpus=20 --dun_use_partial & #30min 300,000 #1h 40w #30cpu 30w 

generate_batch_cmd_for_Cogent_family_finding.py --cpus=20 --cmd_filename=cmd preCluster.cluster_info.csv preCluster_out cDNA

#cogent 分组
bash ~/analysis/scripts/Cogent_check1.sh 


################3.reconstruction#############
printf "Partition\tSize\tMembers\n" > final.partition.txt
ls preCluster_out/*/*partition.txt | xargs -n1 -i sed '1d; $d' {} | cat >> final.partition.txt

#排序 40000的那些bin单独挑出来（如果分了很多个就挑出来）
generate_batch_cmd_for_Cogent_reconstruction.py cDNA/ | sort -n -k1,1 > cmd2





bash ~/analysis/scripts/Cogent_check2.sh |wc -l



################4.fake genome#############
export PATH=$PATH:~/biosoft/cDNA_Cupcake-master/sequence/
tail -n 1 final.partition.txt | tr ',' '\n'  > unassigned.list
get_seqs_from_list.py isoseq_flnc.fasta unassigned.list > unassigned.fasta

mkdir collected
cd collected
cat ../cDNA/*/cogent2.renamed.fasta ../unassigned.fasta > cogent.fake_genome.fasta




################5.gmap & collapse#############
cd collected
# gmap
gmap_build -D . -d fake_genome cogent.fake_genome.fasta #1min 
nohup gmap -D . -d fake_genome -f samse -n 0 -t 30 \
   ../isoseq_flnc.fasta > isoseq_flnc.fasta.sam 2> \
   isoseq_flnc.fasta.sam.log &  #3h30min #30cpu 
#collapse
sort -k 3,3 -k 4,4n isoseq_flnc.fasta.sam > isoseq_flnc.fasta.sorted.sam
nohup collapse_isoforms_by_sam.py --input ../isoseq_flnc.fasta -s isoseq_flnc.fasta.sorted.sam \
           -c 0.95 -i 0.85 --dun-merge-5-shorter \
           -o isoseq_flnc.fasta.no5merge & #5 min

#没有cluster report 以下3步省略（非必须）
#下载cluster report并合并 
ln -s ../../raw/cluster_report.csv
#报错 
get_abundance_post_collapse.py isoseq_flnc.fasta.no5merge.collapsed cluster_report.csv
filter_away_subset.py isoseq_flnc.fasta.no5merge.collapsed


#处理 ID
#去掉"|0_0|path1:1-7569(+)"
#cat cDNA.unique.fa.transdecoder.cds|grep "\|[0-9]+_[0-9]+\|path[0-9]+:[0-9]+-[0-9]+\(.\)" -E -o
cat isoseq_flnc.fasta.no5merge.collapsed.rep.fa|sed 's/|[0-9]\+_[0-9]\+|path[0-9]\+:[0-9]\+-[0-9]\+(.)//g' > cDNA.unique.fa


cd ../../
ln -s cogent/collected/cDNA.unique.fa cDNA.unique.fa


