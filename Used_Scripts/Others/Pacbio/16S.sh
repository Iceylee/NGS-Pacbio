
#download 16S_chemerin_tutorial.zip
#conda安装pear 
# https://github.com/LangilleLab/microbiome_helper
#conda安装fastq-join
  
conda install bbmap vsearch sortmerna picrust #1用于过滤 2去除嵌合体 3去除rna？
MICROPATH=/home/liyubing/bin/microbiome_helper-master

# 1.1合并双端方法1
# join_paired_ends.py -f 105CHE6WT_S325_L001_R1_001.fastq -r 105CHE6WT_S325_L001_R2_001.fastq -m fastq-join -o temp/fastq_join

# for i in `ls *R1_001.fastq`
# do
# 	i=${i/R1_001.fastq/}
# 	echo "join_paired_ends.py -f ${i}R1_001.fastq -r ${i}R2_001.fastq -m fastq-join -o temp/fastq_join/${i}"
# done

# 1.2合并双端方法2:pear
cp /home/liyubing/bin/pear-0.9.11-linux-x86_64/bin/pear /usr/bin/pear
perl run_pear.pl -p 1 -o stitched_reads fastq/*fastq #每对得到4个文件，用assembled.fastq

# 2.fastqc
mkdir fastqc_out_combined
cat stitched_reads/*.assembled.fastq | fastqc -t 1 stdin -o fastqc_out_combined

cd fastqc_out_combined
mv stdin_fastqc.html combined_fastqc.html  
mv stdin_fastqc.zip combined_fastqc.zip  

# 3.过滤:fastx-toolkit bbmap
#修改read_filter.pl 的$bbmap_dir路径/data/software/miniconda2/bin
perl $MICROPATH/read_filter.pl -q 30 -p 90 -l 400 -thread 1 -c both stitched_reads/*.assembled*fastq #输出filter_reads文件夹

# 4.fastqc
mkdir fastqc_out_combined_filtered
cat filtered_reads/*.fastq | fastqc -t 1 stdin -o fastqc_out_combined_filtered

cd fastqc_out_combined_filtered
mv stdin_fastqc.html combined_filtered_fastqc.html
mv stdin_fastqc.zip combined_filtered_fastqc.zip

# 5.转fastq并去除嵌合体
perl $MICROPATH/run_fastq_to_fasta.pl -p 1 -o fasta_files filtered_reads/*fastq
perl $MICROPATH/chimera_filter.pl -type 1 -thread 5 -db $MICROPATH/RDP_trainset16_022016.fa fasta_files/*fasta #non_chimeras
#db 下载This DB was originally from the Ribosome Database Project (RDP) and was parsed to include only bacteria.

# 6.Run open-reference OTU picking pipeline
# 合并所有fasta
add_qiime_labels.py -i non_chimeras/ -m map.txt -c FileInput -o combined_fasta

echo "pick_otus:threads 5" >> clustering_params.txt
echo "pick_otus:sortmerna_coverage 0.8" >> clustering_params.txt
#echo "pick_otus:sortmerna_db /home/shared/pick_otu_indexdb_rna/97_otus" >> clustering_params.txt
pick_open_reference_otus.py -i combined_fasta/combined_seqs.fna -o clustering/ -p clustering_params.txt -m sortmerna_sumaclust -s 0.1 -v --min_otu_size 1   #调用sortmerna
#SortMeRNA for reference picking and SUMACLUST for de novo OTU picking (~24 hours):


# 7.Remove low confidence OTUs
remove_low_confidence_otus.py -i $PWD/clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o $PWD/clustering/otu_table_high_conf.biom
biom summarize-table -i clustering/otu_table_mc1_w_tax_no_pynast_failures.biom -o clustering/otu_table_mc1_w_tax_no_pynast_failures_summary.txt

biom summarize-table -i clustering/otu_table_high_conf.biom -o clustering/otu_table_high_conf_summary.txt


# 8.Rarify reads
mkdir final_otu_tables
single_rarefaction.py -i clustering/otu_table_high_conf.biom -o final_otu_tables/otu_table.biom -d 375

# 9.Diversity analyses











source activate qiime1


export PATH=$PATH:~/biosoft/cDNA_Cupcake-master/sequence/
export PATH=/home/liyubing/miniconda3/bin/perl:$PATH



