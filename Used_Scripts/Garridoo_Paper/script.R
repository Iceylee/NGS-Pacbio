
# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list=ls())

# load paths to project directories
source("paths.R")

# create folder to output figures
system(paste("mkdir -p ", figures.dir, sep=""))

### functional PCoA
# generates Fig. 3A (func. PCoA)
# Fig. 3B (pairwise distances within families) and

source("functional_MDS.R")


### cluster analysis:
# calculate clusters based on AMPHORA genes
# generate supplemental cluster analysis PDF

source("seq_id_clustering.R")

### KEGG pathway analysis
# generages Fig. 4A (func. heatmap & tree of cluster reps.)
# Fig. 4B (enrichment of functions root v. leaf & v. soil) and
# Fig. 4C (pangenome distribution barplots per cluster)

source("func_heatmap_cluster.R")

### tree figures
# generages Fig. S7 (spp. tree + secretion systems)
# generages Fig. S8 (spp. tree + Carbohydrate met. abundances)
# generages Fig. S9 (spp. tree + Xenobiotic met. abundances)

source("tree_figures.R")

