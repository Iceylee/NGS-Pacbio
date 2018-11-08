######################################################################
####### 1.  Get All Seq without chimeric #######
for i in $(find /Databackup/2018_04/ChenXiaoXuan/01.CleanData/ -maxdepth 2  -name "W*.fna")
do echo "usearch7 -uchime_ref $i -db /data/16S/gold.fa -strand plus -nonchimeras `echo $i| awk -F'S/' '{print $2}'|awk -F'.fna' '{print $1.".good.fna"}'`"
done

cat *.good.fna >> AllSeqs.fna 

### 476,922 --> 474,083
######################################################################



######################################################################
####### 2. Pick OTU ########
pick_open_reference_otus.py -i ../1.DeChimeric/AllSeqs.fna -r /data/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta -o otu_output -s 0.1 -m usearch61 -p params.txt

biom convert -i otu_output/otu_table_mc2_w_tax.biom -o otu_table.txt --to-tsv --header-key taxonomy

cp otu_output/otu_table_mc2_w_tax.biom ./otu_table.biom

biom summarize-table -i otu_table.biom -o otu_table_summary.txt
### OTU 978 
### min: 70725, max: 83193 

summarize_taxa_through_plots.py -i otu_table.biom -m ../Metadata.txt -o taxa
######################################################################



######################################################################
######### 3. Get high abundance OTU , and make phylogeny ###########
filter_otus_from_otu_table.py -i otu_table.biom -o otu_table_high_abundance.biom --min_count_fraction 0.00005
##### OTU: 197
##### min: 69979, max: 82774

summarize_taxa_through_plots.py -i otu_table_high_abundance.biom -m ../Metadata.txt -o taxa_high_abundance

biom convert -i otu_table_high_abundance.biom -o otu_table_high_abundance.txt --to-tsv --header-key taxonomy

python /data1/script/ExtraGeneID2Fa.py  tmp.id.txt  otu_output/rep_set.fna  > rep_set_high_abundance.fna 

align_seqs.py -i rep_set_high_abundance.fna -t /data1/software/miniconda2/envs/qiime1/lib/python2.7/site-packages/qiime_test_data/align_seqs/core_set_aligned.fasta.imputed -o pynast_high_abundance_aligned

filter_alignment.py -i pynast_high_abundance_aligned/rep_set_high_abundance_aligned.fasta -o rep_set_high_abundance_filtered_alignment
### 175.fna

make_phylogeny.py -i rep_set_high_abundance_filtered_alignment/rep_set_high_abundance_aligned_pfiltered.fasta -o rep_set_high_abundance.tre
######################################################################




######################################################################
#### 4. HeatMap with OTUs abundance >=1% ####
# ?? # filter_otus_from_otu_table.py -i otu_table_high_abundance.biom -o otu_table_Top30_heatmap.biom --min_count_fraction 0.01 

# ?? # biom summarize-table -i otu_table_Top30_heatmap.biom -o otu_table_Top30_heatmap_summary.txt
## abundance 1% ---> 18 otu  ### 

make_otu_heatmap.py -i otu_table_high_abundance.biom -o otu_table_high_abundance_heatmap.pdf -m ../Metadata.txt 
######################################################################


################### 5. Rarefaction with alpha diversity, and data should enought , This is Just Abundance OTU (197)### 
#重抽样标准化 （single or multiple）
multiple_rarefactions.py -i ../2.PickOTUs/otu_table_high_abundance.biom -m 100 -x 80000 -s 2000 -n 10 -o rarefied_otu_tables/ 

alpha_diversity.py -i rarefied_otu_tables/ -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1  -o rarefied_alpha_diversity/ -t ../2.PickOTUs/rep_set_high_abundance.tre

collate_alpha.py -i rarefied_alpha_diversity -o alpha_collated

make_rarefaction_plots.py -i alpha_collated/ -m ../Metadata.txt -o rarefaction_plots
######################################################################


####################### Alpha diversity with all OTU #####################
multiple_rarefactions.py -i ../2.PickOTUs/otu_table.biom -m 100 -x 80000 -s 2000 -n 10 -o rarefied_otu_tables

alpha_diversity.py -i rarefied_otu_tables/ -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o rarefied_AlphaDiversity_all/ -t ../2.PickOTUs/otu_output/rep_set.tre

collate_alpha.py -i rarefied_AlphaDiversity_all/ -o alpha_collated_all

make_rarefaction_plots.py -i alpha_collated_all/ -m ../Metadata.txt  -o rarefied_plots

#############################################################################

alpha_diversity.py -i ../2.PickOTUs/otu_table.biom -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o alpha_diversity -t ../2.PickOTUs/otu_output/rep_set.tre 


######################################################################
#### 6. Alpha diversity ###
alpha_diversity.py -i ../2.PickOTUs/otu_table_high_abundance.biom -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o rarefied_alpha_diversity -t ../2.PickOTUs/rep_set_high_abundance.tre

######################################################################

#### mv rarefied_alpha_diversity 4.Alpha/rarefied_alpha_diversity_otu_table
#### collate_alpha.py -i ../4.Alpha/rarefied_alpha_diversity_otu_table -o alpha_collated_otu_table
####

######################################################################
#### 7. Beta diversity  and PCoA plots ###  Should normalized ####
normalize_table.py -i ../2.PickOTUs/otu_table.biom  -a CSS -o ../5.Beta/normalized_otu_table.biom

beta_diversity.py -i normalized_otu_table.biom -m weighted_unifrac,unweighted_unifrac -o beta_div_normalized -t ../2.PickOTUs/otu_output/rep_set.tre

principal_coordinates.py -i beta_div_normalized/ -o pcoa_normalized

make_2d_plots.py -i pcoa_normalized/ -m ../Metadata.txt -o 2d_plots_normalized


######################################################################
#### 7. 2d PCoA, Treatment Distance, ###
make_distance_boxplots.py -d ../5.Beta/beta_div/weighted_unifrac_otu_table_high_abundance.txt -m ../Metadata.txt -o Treatment_distance -f 'Treatment'

compare_categories.py --method anosim -i ../5.Beta/beta_div/weighted_unifrac_otu_table_high_abundance.txt -m ../Metadata.txt  -c 'Treatment' -n 1000 -o anosim_out 

######################################################################



differential_abundance.py -i ../2.PickOTUs/otu_table_high_abundance.biom  -o DiffAbundance_high_otus.txt -m ../Metadata.txt -a DESeq2_nbinom -c Treatment -x Treat -y CK -d 



differential_abundance.py -i ../2.PickOTUs/otu_table.biom  -o DiffAbundance_all_otus.txt -m ../Metadata.txt -a DESeq2_nbinom -c Treatment -x Treat -y CK -d

###
### 3D PCoA plots with Emperor ### 
jackknifed_beta_diversity.py -i ../2.PickOTUs/otu_table.biom -o bdiv_jk100 -e 100 -m ../Metadata.txt -t ../2.PickOTUs/otu_output/rep_set.tre 

dissimilarity_mtx_stats.py -i bdiv_jk100/unweighted_unifrac/rare_dm/ -o dissimilar_stat_output 


### between and within Treatment's distance ###
make_distance_boxplots.py -m ../Metadata.txt -o boxplot_output -d dissimilar_stat_output/means.txt -f 'Treatment' --save_raw_data


#### Show similarity of bacterial communities based on 16s rRNA genes ###  Weight & Unweight  ######
make_bootstrapped_tree.py -m bdiv_jk100/unweighted_unifrac/upgma_cmp/master_tree.tre -s bdiv_jk100/unweighted_unifrac/upgma_cmp/jackknife_support.txt -o bdiv_jk100/unweighted_unifrac/upgma_cmp/Tree.pdf

make_bootstrapped_tree.py -m bdiv_jk100/weighted_unifrac/upgma_cmp/master_tree.tre -s bdiv_jk100/weighted_unifrac/upgma_cmp/jackknife_support.txt -o bdiv_jk100/weighted_unifrac/upgma_cmp/Tree.pdf


####### R: PCA  otu_table_high_abundance.txt ################ 
abundance_data.pca <- princomp(abundance_data[,1:6])
plot3d(abundance_data.pca$loadings[,1:3], col=colour, type = "s", radius=0.025) + legend3d("topright", legend=colnames(abundance_data[,1:6]), col=colour, pch=16)
