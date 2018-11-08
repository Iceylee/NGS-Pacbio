source activate qiime1

# 1.DropChimeric

#usearch7 -uchime_ref /Databackup/2018_04/ChenXiaoXuan/01.CleanData/WR180386S/WR180386S.fna -db /data/16S/gold.fa -strand plus -nonchimeras WR180386S.good.fna

for smp in WR180384S WR180385S WR180386S WR180387S WR180388S WR180389S 
do
usearch7 -uchime_ref clean_data/$smp.fna -db /data1/16S/gold.fa -strand plus -nonchimeras 1.Assemble/$smp.good.fna
done


cat *.good.fna >> AllSeqs.fna

######################################################################
####### 2. Pick OTU ########
pick_open_reference_otus.py -i 1.Assemble/AllSeqs.fna -r /data1/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -o 2.PickOTUs/otu_output -s 0.1 -m usearch61 -p Params.txt

cd 2.PickOTUs
biom convert -i otu_output/otu_table_mc2_w_tax.biom -o otu_table.txt --to-tsv --header-key taxonomy

cp otu_output/otu_table_mc2_w_tax.biom ./otu_table.biom

biom summarize-table -i otu_table.biom -o otu_table_summary.txt

summarize_taxa_through_plots.py -i otu_table.biom -m ../MetaData.txt -o taxa

######################################################################
######### 3. 获取高丰度OTU , and make phylogeny系统进化树 ###########
filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_high_abundance.biom --min_count_fraction 0.00005
##### OTU: 197
##### min: 69979, max: 82774

summarize_taxa_through_plots.py -i otu_table_high_abundance.biom -m ../MetaData.txt -o taxa_high_abundance

biom convert -i otu_table_high_abundance.biom -o otu_table_high_abundance.txt --to-tsv --header-key taxonomy

python /data1/script/Pipeline/16S/ExtraOTU2Seq.py otu_table_high_abundance.txt otu_output/rep_set.fna > rep_set_high_abundance.fna 


align_seqs.py -i rep_set_high_abundance.fna -t /data1/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_test_data/align_seqs/core_set_aligned.fasta.imputed -o pynast_high_abundance_aligned

#筛选结果中保守序列和保守区
filter_alignment.py -i pynast_high_abundance_aligned/rep_set_high_abundance_aligned.fasta -o rep_set_high_abundance_filtered_alignment

#构建进化树
make_phylogeny.py -i rep_set_high_abundance_filtered_alignment/rep_set_high_abundance_aligned_pfiltered.fasta -o rep_set_high_abundance.tre

################################
#OTU 热图 heatmap
################################
make_otu_heatmap.py -i otu_table_high_abundance.biom -o otu_table_high_abundance.pdf -m ../MetaData.txt


####################### 
#Alpha diversity with all OTU 
#####################
#重抽样标准化 （single or multiple）
multiple_rarefactions.py -i ../2.PickOTUs/otu_table.biom -m 100 -x 80000 -s 2000 -n 10 -o rarefied_otu_tables

# 计算常用的四种Alpha多样性指数
alpha_diversity.py -i rarefied_otu_tables/ \
	-m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 \
	-o rarefied_AlphaDiversity_all/ \
	-t ../2.PickOTUs/otu_output/rep_set.tre 

#得到各参数数据
collate_alpha.py -i rarefied_AlphaDiversity_all/ -o alpha_collated_all
#作图
make_rarefaction_plots.py -i alpha_collated_all/ -m ../MetaData.txt  -o rarefied_plots

#############################################################################
#不做重抽样标准化 -i
alpha_diversity.py -i ../2.PickOTUs/otu_table.biom -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o alpha_diversity -t ../2.PickOTUs/otu_output/rep_set.tre 

alpha_diversity.py -i ../2.PickOTUs/otu_table_high_abundance.biom -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o rarefied_alpha_diversity -t ../2.PickOTUs/rep_set_high_abundance.tre




################################
#### 7. Beta diversity  and PCoA plots 
###  Should normalized ####
######################################
cd 6.Beta_diversity
#CSS 标准化
normalize_table.py -i ../2.PickOTUs/otu_table.biom  -a CSS -o normalized_otu_table.biom

beta_diversity.py -i normalized_otu_table.biom -m weighted_unifrac,unweighted_unifrac -o beta_div_normalized -t ../2.PickOTUs/otu_output/rep_set.tre

principal_coordinates.py -i beta_div_normalized/ -o pcoa_normalized

make_2d_plots.py -i pcoa_normalized/ -m ../MetaData.txt -o 2d_plots_normalized


###################################
#### 7. 2d PCoA, Treatment Distance, 
###################################
make_distance_boxplots.py -d beta_div_normalized/weighted_unifrac_otu_table_high_abundance.txt -m ../MetaData.txt -o Treatment_distance -f 'Treatment'

compare_categories.py --method anosim -i ../5.Beta/beta_div/weighted_unifrac_otu_table_high_abundance.txt -m ../Metadata.txt  -c 'Treatment' -n 1000 -o anosim_out 

