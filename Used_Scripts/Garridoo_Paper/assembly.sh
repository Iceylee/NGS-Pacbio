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

if [ ! -f $config_file ]
then
echo "invalid config file"
    return 1
fi

# load config file
source $config_file

# load functions
source assembly.functions.sh

for genome_id in $(cut -f 2 $mapping_file | tail -n +2)
do

    log "["$genome_id"] processing raw data..."

    zcat $data_dir/"$genome_id"_1.fastq.gz >> $working_dir/"$genome_id"_1.fastq
    zcat $data_dir/"$genome_id"_2.fastq.gz >> $working_dir/"$genome_id"_2.fastq

    trim -genome_id=$genome_id -mode=$mode -lead=$lead \
         -trail=$trail -sw_l=$sw_length -sw_q=$sw_min_q \
         -min_l=$min_length -n_cores=$n_cores \
         &>> $output

    assemblySOAP -genome_id=$genome_id -mode=$mode \
                 -kmer=$kmer -avg_ins=$avg_ins \
                 -max_rd_len=$max_rd_len -n_cores=$n_cores \
                 &>> $output

    assemblyA5 -genome_id=$genome_id -n_cores=$n_cores \
                &>> $output

    assemblyStats $genome_id

    # cleanup

    rm -f $working_dir/"$genome_id"_*fastq
    rm -f -r $working_dir/SOAP/"$genome_id"
    rm -f -r $working_dir/A5/"$genome_id"

    log "["$genome_id"] done"

done

log "DONE!"

