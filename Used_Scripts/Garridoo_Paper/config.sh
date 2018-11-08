#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# paths

# output and intermediate data folder
working_dir="/net/phylogenetics/atsphere/working_dir"
# intermediate data folder
data_dir="/net/phylogenetics/atsphere/data"
# reference databases folder
path_to_db="/net/phylogenetics/atsphere/refdata"
# mapping file with isolates medatada
mapping_file="$data_dir/mapping.txt"

# debug
logfile="/dev/tty"
output=""$working_dir"/output.txt"

# parameters (general)
n_cores=$(nproc)

# parameters (trimmomatic)
min_length=36
lead=3
trail=3
sw_length=4
sw_min_q=15

# parameters (SOAPdenovo)
assembler="SOAP"
mode="PE"
kmer=63
avg_ins=350
max_rd_len=100

# parameters (MrBayes)
ngen=30000

# paremeters (annotation)
k=9
rel=2
E="1.0e-5"
min_cov=70

# parameters (taxator)
taxator_speedup=0.9
refpack="/tmp/nonredundant-microbial"

