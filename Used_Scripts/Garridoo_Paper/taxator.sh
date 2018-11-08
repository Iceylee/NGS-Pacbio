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

genomes=$(cut -f 1 $data_dir/assembly_table.csv | tail -n +2)

mkdir -p $working_dir/taxator/logs
cd $working_dir/taxator

# serial

for i in $genomes
do

    echo $(date -u) $i
    rm -f $working_dir/taxatir/logs/"$i".log

    binning-last.bash $refpack $data_dir/assemblies/"$i".fna &>> $working_dir/taxator/logs/"$i".log

done

log "DONE!"

