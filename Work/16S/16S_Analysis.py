#!/usr/bin/env python
# coding: utf-8

### Author: zss
### Date: 2018.7.16
## For: 16S Analysis ###

import ConfigParser
import sys
import os
import commands

if len(sys.argv) != 2:
    print "\n\t\tUsage: python %s 16S_Analysis.conf\n" % sys.argv[0]
    sys.exit(1)

ConfgFile = sys.argv[1]
Config = ConfigParser.ConfigParser()
Config.read(ConfgFile)

### DB
gold_16S_Seq = Config.get("DB", "gold_16S_Seq")
gg_13_8_otus = Config.get("DB", "gg_13_8_otus")

### Get Data Info
DataDir = Config.get("Data", "Dir")
OutPut = Config.get("Data", "OutPut")

### Get Qiime1
MetaData = Config.get("Qiime1", "MetaData")
ParamsFile = Config.get("Qiime1", "ParamsFile")
min_count_fraction = Config.get("Qiime1", "min_count_fraction")

#####################################################################################################################
print "\n################################################ Info ###############################################################\n"
print "\t\tData Dir: %s"  % DataDir
print "\t\tOutPut Dir: %s"  % OutPut 
print "\t\tmin_count_fraction: %s"  % min_count_fraction
print "\n#####################################################################################################################\n"
#####################################################################################################################


def MkdirOutPut(OutPut):
    for i in [ "0.QC", "1.Assemble", "2.PickOTUs", "3.Taxonomy", "4.Rarefaction", "5.Alpha_diversity", "6.Beta_diversity", "7.PCA_PCoA/pcoa_normalized", "7.PCA_PCoA/2d_plots_normalized", "8.Other_Analysis"]:
        if not os.path.exists(OutPut + i):
            os.makedirs(OutPut + i)


def DropChimeric(DataDir, OutPut, gold_16S_Seq):
    Command1 = "for i in $(find %s -maxdepth 2  -name 'W*.fna');do usearch7 -uchime_ref $i -db %s -strand plus -nonchimeras `echo $i| awk -F'S/' '{print $2}'| awk -F'.fna' '{print $1.'.good.fna'}'`;done" % (DataDir, gold_16S_Seq)

    print "\033[1;31;40m --** Filter Data **--Command1: \033[0m \n%s" % (Command1)
    os.system(Command1)
            
    Command2 = "cat *.good.fna >> AllSeqs.fna"
    print "\033[1;34;40m --** Cat all good fna to AllSeqs.fna **--Command2: \033[0m \n%s" % Command2
    os.system(Command2)


def PickOTU(AllSeqFile, MetaData, gg_13_8_otus, OutPut, ParamsFile):
    OTU_Dir = OutPut + "2.PickOTUs/otu_output"

    ### pick otu ###
    Command3 = "pick_open_reference_otus.py -i %s -r %s -o %s -p %s \n\n" % (AllSeqFile, gg_13_8_otus, OTU_Dir, ParamsFile)
    print "--** OTU Pick **--Command3: \n%s" % Command3
    os.system(Command3)

    OTU_Table_biom = OTU_Dir + "otu_table_mc2_w_tax_no_pynast_failures.biom"
    OTU_Table_biom_alias = OTU_Dir + "otu_table.biom"
    OTU_Table_txt = OTU_Dir + "otu_table.txt"

    ## alias ###
    Command3_1 = "ln -S %s %s" % (OTU_Table_biom, OTU_Table_biom_alias)
    os.system(Command3_1)

    ### otu biom to txt ###
    Command4 = "biom convert -i %s -o %s --to-tsv --header-key taxonomy" % (OTU_Table_biom, OTU_Table_txt)
    print "--** OTU biom to txt **--Command4: \n%s\n\n" % Command4
    os.system(Command4)

    ### otu biom to summary ###
    OTU_Table_summary = OTU_Dir + "otu_table_summary.txt"
    Command5 = "biom summarize-table -i %s -o %s" % (OTU_Table_biom, OTU_Table_summary)
    print "--** OTU biom to summary **--Command5: \n%s\n\n" % Command5
    os.system(Command5)

    ### otu taxa and plot ###
    Taxa_Dir = OutPut + "3.Taxonomy/" 
    Command6_1 = "summarize_taxa_through_plots.py -i %s -m %s -o %s" % (OTU_Table_biom, MetaData, Taxa_Dir)
    print "--** OTU summarize taxa and plot **--Command6_1: \n%s\n\n" % Command6_1
    os.system(Command6_1)

    ### otu network ###
    OTU_Network = OutPut + "3.Taxonomy/OTU_Network"
    Command6_2 = "make_otu_network.py -i %s -m %s -o %s" % (OTU_Table_biom, MetaData, OTU_Network)
    print "--** OTU Network **--Command6_2: \n%s\n\n" % Command6_2
    os.system(Command6_2)


def Get_High_Abundance_OTU(OutPut, min_count_fraction, MetaData):
    OTU_Table_biom_alias = OutPut + "2.PickOTUs/otu_output/otu_table.biom"
    High_OTU_Table_biom = OutPut + "2.PickOTUs/HighAbundance/otu_table_high_abundance.biom"
    High_OTU_Table_txt = OutPut + "2.PickOTUs/HighAbundance/otu_table_high_abundance.txt"

    rep_set_fa = OutPut + "2.PickOTUs/otu_output/rep_set.fna"
    High_rep_set_fa = OutPut + "2.PickOTUs/HighAbundance/rep_set_high_abundance.fna"

    Command7 = "filter_otus_from_otu_table.py -i %s -o %s --min_count_fraction %s\n\n" % (OTU_Table_biom_alias, High_OTU_Table_biom, min_count_fraction)
    print "--** Get High Abundance OTU **--Command7: \n%s" % Command7
    os.system(Command7)

    ### high otu biom to txt ###  
    Command8 = "biom convert -i %s -o %s --to-tsv --header-key taxonomy" % (High_OTU_Table_biom, High_OTU_Table_txt)
    print "--** OTU biom to txt **--Command8: \n%s" % Command8
    os.system(Command8)

    ### otu taxa and plot ###
    High_Taxa_Dir = OutPut + "2.PickOTUs/HighAbundance/Taxonomy/"
    Command9 = "summarize_taxa_through_plots.py -i %s -m %s -o %s" % (OTU_Table_biom_HighAbundance, MetaData, High_Taxa_Dir)
    print "--** High OTU summarize taxa and plot **--Command9: \n%s" % Command9
    os.system(Command9)

    ### get high abundance seq ### 
    Command10 = "python /data1/script/Pipeline/16S/ExtraOTU2Seq.py %s %s > %s" % (High_OTU_Table_txt,rep_set_fa, High_rep_set_fa )
    print "--** High OTU seq **--Command10: \n%s\n\n" % Command10
    os.system(Command10)
    
    Command11 = "make_otu_heatmap.py -i %s -o %s -m %s" % (High_OTU_Table_biom, High_OTU_Table_HeatMap, MetaData)
    print "--** High OTU HeatMap **--Command11: \n%s\n\n" % Command11
    os.system(Command11)


def Alpha_Diversity(OutPut, MetaData):
    OTU_Table_biom = OutPut + "2.PickOTUs/otu_output/otu_table.biom"
    rep_set_tre = OutPut + "2.PickOTUs/otu_output/rep_set.tre"
    alpha_diversity_Dir = OutPut + "5.Alpha_diversity/"

    Command13_1 = "alpha_diversity.py -i %s -m observed_species,shannon,PD_whole_tree,singles,simpson,observed_otus,chao1 -o %s -t %s" % (OTU_Table_biom, alpha_diversity_Dir, rep_set_tre)
    print "--** aplha diversity **--Command13_1: \n%s\n\n" % Command13_1
    os.system(Command13_1)

    #Command13_2 = "alpha_rarefaction.py -i %s -m %s -p %s -t %s -e %s" % (OTU_Table_biom, alpha_diversity_Dir, rep_set_tre)
    #print "--** aplha rarefaction **--Command13_2: \n%s\n\n" % Command13_2
    #os.system(Command13_2)


def Beta_Diversity(OutPut, MetaData):
    OTU_Table_biom = OutPut + "2.PickOTUs/otu_output/otu_table.biom"
    rep_set_tre = OutPut + "2.PickOTUs/otu_output/rep_set.tre"
    Normalized_OTU_Table_biom = OutPut + "6.Beta_diversity/normalized_otu_table.biom"
    
    Command14 = "normalize_table.py -i %s -a CSS -o %s" % (OTU_Table_biom, Normalized_OTU_Table_biom)
    print "--** beta diversity: normalize otu **--Command14: \n%s\n\n" % Command14
    os.system(Command14)
 
    beta_diversity_Dir = OutPut + "6.Beta_diversity/"
    Command15 = "beta_diversity.py -i %s -m weighted_unifrac,unweighted_unifrac -o %s -t %s" % (Normalized_OTU_Table_biom, beta_diversity_Dir, rep_set_tre)
    print "--** beta diversity **--Command15: \n%s\n\n" % Command15
    os.system(Command15)

    ### PCOA ###
    pcoa_Dir = OutPut + "7.PCA_PCoA/pcoa_normalized/"
    Command16 = "principal_coordinates.py -i %s -o %s" % (beta_diversity_Dir, pcoa_Dir)
    print "--** beta diversity: PCOA **--Command16: \n%s\n\n" % Command16
    os.system(Command16)

    ### 2d plot ### 
    pcoa_2d_plot_Dir = OutPut + "7.PCA_PCoA/2d_plots_normalized/"
    Command17 = "make_2d_plots.py -i %s -m %s -o %s" % (pcoa_Dir, MetaData, pcoa_2d_plot_Dir)
    print "--** beta diversity: 2d plot **--Command17: \n%s\n\n" % Command17
    os.system(Command17)

    ### 3d PCoA ###
    bdiv_jk100 = OutPut + "8.Analysis/bdiv_jk100/"
    Command18 = "jackknifed_beta_diversity.py -i %s -o %s -e 100 -m %s -t %s" % (OTU_Table_biom, bdiv_jk100, MetaData, rep_set_tre)
    print "--** jackknifed_beta_diversity **--Command18: \n%s\n\n" % Command18
    os.system(Command18)

    ### Show similarity of bacterial communities based on 16s rRNA genes ###  Weight & Unweight ###
    ### Unweight ###
    bdiv_jk100_master_tree = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/master_tree.tre"
    jackknife_support = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/jackknife_support.txt"
    unweight_Tree = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/Tree.pdf"
    Command19 = "make_bootstrapped_tree.py -m %s -s %s -o %s" % (bdiv_jk100_master_tree, jackknife_support, unweight_Tree)

    ### Weight ###
    wei_bdiv_jk100_master_tree = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/master_tree.tre"
    wei_jackknife_support = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/jackknife_support.txt"
    weight_Tree = OutPut + "8.Analysis/bdiv_jk100/unweighted_unifrac/upgma_cmp/Tree.pdf"
    Command20 = "make_bootstrapped_tree.py -m %s -s %s -o %s" % (wei_bdiv_jk100_master_tree, wei_jackknife_support, weight_Tree)
    print "--** Unweight and Weight Tree  **--Command19: %s\n\n Command20: \n%s\n\n" % (Command19, Command20)
    os.system(Command19)
    os.system(Command20)


def 


if __name__ == '__main__':
    print "\t\t\033[1;32;40m ### Step1: Mkdir ###\033[0m \n\n"
    MkdirOutPut(OutPut)

    print "\t\t\033[1;32;40m ### Step2: Drop Chimeric Seq ###\033[0m \n\n"
    DropChimeric(DataDir, OutPut, gold_16S_Seq)
 
    print "\t\t\033[1;32;40m ### Step3: Pick OTU ###\033[0m \n\n"
    PickOTU(AllSeqFile, MetaData, gg_13_8_otus, OutPut, ParamsFile)
 
    print "\t\t\033[1;32;40m ### Step4: Get High Abundance OTU ###\033[0m \n\n"
    Get_High_Abundance_OTU(OutPut, min_count_fraction, MetaData)

    print "\t\t\033[1;32;40m ### Step5: Alpha_Diversity ###\033[0m \n\n"
    Alpha_Diversity(OutPut)

    print "\t\t\033[1;32;40m ### Step6: Beta_Diversity ###\033[0m \n\n"
    Beta_Diversity(OutPut, MetaData)
