#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

log() {

    echo $(date -u)": "$1 >> $logfile

}

config () {

    # parse arguments
    local config_file=$1

	# check config file
	if [ ! -f $config_file ]
	then
	    echo "invalid config file"
	    return 1
	fi
	
	# get dir and file paths
	source $config_file

}

trim () {

    # default parameters
    local n_cores=$(nproc)
    local lead=3
    local trail=3
    local sw_l=4
    local sw_q=15
    local min_l=36
    local mode="PE"

    # parse arguments
    for i in "$@"
    do
        case $i in
            -n_cores=* )
            n_cores=`echo $i | sed 's/^.n_cores.//g'` ;;
            -genome_id=* )
            genome_id=`echo $i | sed 's/^.genome_id.//g'` ;;
            -lead=* )
            lead=`echo $i | sed 's/^.lead.//g'` ;;
            -trail=* )
            trail=`echo $i | sed 's/^.trail.//g'` ;;
            -sw_l=* )
            sw_l=`echo $i | sed 's/^.sw_l.//g'` ;;
            -sw_q=* )
            sw_q=`echo $i | sed 's/^.sw_q.//g'` ;;
            -min_l=* )
            min_l=`echo $i | sed 's/^.min_l.//g'` ;;
            -mode=* )
            mode=`echo $i | sed 's/^.mode.//g'` ;;
        esac
    done
    
    if [[ (! -f $data_dir/"$genome_id"_1.fastq || ! -f $data_dir/"$genome_id"_2.fastq) &&  ! -f $data_dir/"$genome_id".fastq ]]
    then
        log "[genome "$genome_id"] fastq file not found!"
        return 1
    fi

    rm -f $working_dir/"$genome_id"_forward_paired.fastq $working_dir/"$genome_id"_1_forward_unpaired.fastq \
    	  $working_dir/"$genome_id"_reverse_paired.fastq $working_dir/"$genome_id"_2_reverse_unpaired.fastq

    log "[genome "$genome_id"] trimming raw reads..."

    if [ "$mode" = "PE" ]
    then
    
        java -jar $trim_path/trimmomatic-0.32.jar PE \
    	     -threads $n_cores \
             -phred33 $data_dir/"$genome_id"_1.fastq $data_dir/"$genome_id"_2.fastq \
    	     $working_dir/"$genome_id"_forward_paired.fastq $working_dir/"$genome_id"_1_forward_unpaired.fastq \
    	     $working_dir/"$genome_id"_reverse_paired.fastq $working_dir/"$genome_id"_2_reverse_unpaired.fastq \
    	     MINLEN:$min_length \
    	     ILLUMINACLIP:$trim_path/adapters/TruSeq3-PE-2.fa:2:40:15 \
    	     LEADING:$lead \
    	     TRAILING:$trail \
             MINLEN:$min_length \
             &>> $output
             # SLIDINGWINDOW:$sw_length:$sw_min_q \

    else if [ "$mode" = "SE" ]
    then

         java -jar $trim_path/trimmomatic-0.32.jar SE \
    	     -threads $n_cores \
             -phred33 $data_dir/"$genome_id".fastq \
    	     $working_dir/"$genome_id".fastq \
    	     MINLEN:$min_length \
    	     ILLUMINACLIP:$trim_path/adapters/TruSeq3-SE.fa:2:40:15 \
    	     LEADING:$lead \
    	     TRAILING:$trail \
    	     MINLEN:$min_length \
             &>> $output 
            #SLIDINGWINDOW:$sw_length:$sw_min_q \

    fi
    fi

}

assemblyA5 () {
    
    # default parameters 
    local n_cores=$(nproc)

    # parse arguments
    for i in "$@"
    do
        case $i in
            -n_cores=* )
            n_cores=`echo $i | sed 's/^.n_cores.//g'` ;;
            -genome_id=* )
            genome_id=`echo $i | sed 's/^.genome_id.//g'` ;;
        esac
    done
    
    #~ # prepare folder for assembly file

    log "[genome "$genome_id"] initializing assembly directory..."
    
    local assembly_dir="$working_dir"/A5/"$genome_id"
    if [ -d "$assembly_dir" ]
    then
        local ndatfiles=$(ls -A "$assembly_dir" | grep $genome_id | wc -l)
        if [ $ndatfiles -gt 0 ]
        then
            log "[genome "$genome_id"] assembly folder is not empty!"
            return 1
        fi
    else
        mkdir -p $assembly_dir
    fi
    
    # run the assembler

    log "[genome "$genome_id"] running A5 assembler..."
    
    local dir=$PWD

    cd $assembly_dir

    a5_pipeline.pl --threads=$n_cores \
                   "$working_dir"/"$genome_id"_forward_paired.fastq \
                   "$working_dir"/"$genome_id"_reverse_paired.fastq \
                   $genome_id

   cd $dir
   
}

assemblySOAP () {

    # default parameters 
    local n_cores=$(nproc)
    local max_rd_len=100
    local avg_ins=200
    local kmer=63
    local mode="PE"

    # parse arguments
    for i in "$@"
    do
        case $i in
            -n_cores=* )
            n_cores=`echo $i | sed 's/^.n_cores.//g'` ;;
            -genome_id=* )
            genome_id=`echo $i | sed 's/^.genome_id.//g'` ;;
            -max_rd_len=* )
            max_rd_len=`echo $i | sed 's/^.max_rd_len.//g'` ;;
            -avg_ins=* )
            avg_ins=`echo $i | sed 's/^.avg_ins.//g'` ;;
            -kmer=* )
            kmer=`echo $i | sed 's/^.kmer.//g'` ;;
            -mode=* )
            mode=`echo $i | sed 's/^.mode.//g'` ;;
        esac
    done
    
    # prepare folder for assembly file

    log "[genome "$genome_id"] initializing assembly directory and creating config files..."
    
    local assembly_dir="$working_dir"/SOAP/"$genome_id"
    if [ -d "$assembly_dir" ]
    then
        local ndatfiles=$(ls -A "$assembly_dir" | grep $genome_id | wc -l)
        if [ $ndatfiles -gt 0 ]
        then
            log "[genome "$genome_id"] assembly folder is not empty!"
            return 1
        fi
    else
        mkdir -p $assembly_dir
    fi
    
    # prepare config file

    if [ "$mode" = "PE" ]
    then

        cat soap_template.config > $assembly_dir/"$genome_id".config
        echo "q1="$working_dir"/"$genome_id"_forward_paired.fastq" >> $assembly_dir/"$genome_id".config
        echo "q2="$working_dir"/"$genome_id"_reverse_paired.fastq" >> $assembly_dir/"$genome_id".config
        echo "q="$working_dir"/"$genome_id"_1_forward_unpaired.fastq" >> $assembly_dir/"$genome_id".config
        echo "q="$working_dir"/"$genome_id"_2_reverse_unpaired.fastq" >> $assembly_dir/"$genome_id".config

    else if [ "$mode" = "SE" ]
    then

        cat soap_template.config > $assembly_dir/"$genome_id".config
        echo "q="$working_dir"/"$genome_id".fastq" >> $assembly_dir/"$genome_id".config

    fi
    fi

    sed -i "s/max_rd_len=/&"$max_rd_len"/g" $assembly_dir/"$genome_id".config
    sed -i "s/avg_ins=/&$avg_ins/g" $assembly_dir/"$genome_id".config

    # run the assembler

    log "[genome "$genome_id"] running SOAP assembler..."
    
    SOAPdenovo-63mer all -s $assembly_dir/"$genome_id".config \
                         -K $kmer -R -p $n_cores \
                         -o $assembly_dir/"$genome_id" \
                         &>$assembly_dir/"$genome_id".log
   
}

assemblyStats () {

    local genome_id=$1

    log "[genome "$genome_id"] generating assembly statistics..."

    $PWD/assembly.stats.R $working_dir/SOAP/$genome_id/"$genome_id".scafSeq \
                                       $genome_id \
                                       $working_dir/SOAP \
                                       "SOAP" \
                                       &> $output

    $PWD/assembly.stats.R $working_dir/A5/$genome_id/"$genome_id".final.scaffolds.fasta \
                                       $genome_id \
                                       $working_dir/A5 \
                                       "A5" \
                                       &> $output

}

