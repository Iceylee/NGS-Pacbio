#usearch 下载
http://www.drive5.com/usearch/download.html
chmod +x usearch6.0.98_i86linux32


#16S 
/data2/ClientData/2018_04/ChenXiaoXuan/Qiime1_4.28

0.FastQC                     
1.DeChimeric                
2.PickOTUs
3.Rarefy
4.Alpha
5.Beta    
6.PCoA 
7.Analysis

###

export PATH=$PATH:/biostack/tools/pipelines/div_seq-0.2.3/bin/utils

export PATH=$PATH:/biostack/labs/seqtk_utils/release

export PATH=$PATH:/biostack/labs/tabtk_utils/release
export PATH=$PATH:/biostack/tools/microbiome/KronaTools-2.7/bin

#
conda install seqtk Trimal mafft fasttree 
seqtk_utils ？ fastx-utils
tabtk_utils
biom 
KronaTools 
phylommand

git clone https://github.com/jameslz/fastx-utils.git
git clone https://github.com/mr-y/phylommand.git
#下载https://github.com/marbl/Krona/wiki/Downloads


