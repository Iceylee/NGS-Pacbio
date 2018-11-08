#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# exits whenever a function returns 1
set -e
# exits if unset variables are used
set -o nounset

# parse arguments
config_file=$1

# check paths file
if [ ! -f $config_file ]
then
    echo "invalid config file"
    return 1
fi

# load paths
source $config_file

# load functions
source datagen.functions.sh

genomes=$(cut -f 1 $mapping_file | tail -n +2)

# process genomes (serial)

for genome_id in $genomes

    processGenome -id $genome_id \
                  -f $data_dir/"$genome_id".fna \
                  -o $working_dir \
                  -n $n_cores \
                  -E $E \
                  -min_cov $min_cov \
                  -k $k \
                  -r $rel \
                  -db ""$path_to_db"/kegg" \
                  &>> $output


done

# organize data

for genome_id in $genomes
do

    mkdir -p $working_dir/assemblies \
             $working_dir/ORFs \
             $working_dir/annotations \
             $working_dir/secretomes

    mv $working_dir/"$genome_id".fna $working_dir/assemblies
    mv $working_dir/"$genome_id".faa $working_dir/ORFs
    mv $working_dir/"$genome_id".ffn $working_dir/ORFs
    mv $working_dir/"$genome_id".ko $working_dir/annotations
    mv $working_dir/"$genome_id"_KEGG.txt $working_dir/annotations
    mv $working_dir/"$genome_id"_SEED.txt $working_dir/annotations

done

log "DONE!"

