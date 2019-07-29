#from snakemake.utils import R,validate, min_version
from snakemake.utils import R,min_version
import sys
import os
import pandas as pd
##### set minimum snakemake version #####
#min_version("5.1.2")

"""
Author: Icey
Date: 2018/9/5 
Run: snakemake   -s Snakefile   
Latest modification: 
"""

##-----------------------------------------------##
## load config and sample sheets 
##-----------------------------------------------##

configfile: "config.yaml"


##-----------------------------------------------##
## A set of functions
##-----------------------------------------------##

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")

##-----------------------------------------------##
## Working directory                             ##
##-----------------------------------------------##

message("The current working directory is " + config["workdir"])


##--------------------------------------------------------------------------------------##
## target rules
##--------------------------------------------------------------------------------------##

rule all:
    input:"outData/prodigal/genome.pep","outData/diamond/sp.anno.out","outData/plots/COG_barplot.pdf","outData/plots/GO_barplot.pdf","outData/diamond/cog_count_fun.txt","outData/diamond/kegg.anno.out","outData/plots/KEGG_barplot.pdf","outData/repeatmasker/RepeatMasker.txt","outData/trf/Trf.txt","outData/sRNA/cmscan.out","outData/sRNA/sRNA.gff","outData/tRNA/tRNA.gff","outData/rRNA/rRNA.gff","outData/island/all.GIs.out","outData/island/GeneIslands.txt","outData/transposon/Transposon.txt","outData/crispr2/Crisprs.txt","outData/anno/Annotation_Summary.txt","outData/prophage/prophage.txt","outData/circos/circos.png"
##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####
include: "rules/1.CDS_Function.smk"
include: "rules/2.Structure.smk"
include: "rules/3.circos.smk"
