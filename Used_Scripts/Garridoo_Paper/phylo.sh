#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e

# parse arguments
config_file="$1"

# check paths file
if [ ! -f $config_file ]
then
    echo "invalid config file"
    return 1
fi

# load paths
source $config_file

# load functions
source phylo.functions.sh

# cleanup
rm -f $working_dir/amphora_all.ffn

# extract sequences, align and build tree

genome_list=$(cut -f 1 $data_dir/mapping.txt | grep -v "ID")
amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

# for each genome in the mapping file

for genome_id in $genome_list 
do

    getAmphoraGenes $genome_id

done

alignAmphoraSeqs
concatenateAmphoraSeqs
buildMLTree

./phylo.R $mapping_file $working_dir/amphora.tree $working_dir 

log "DONE!"

