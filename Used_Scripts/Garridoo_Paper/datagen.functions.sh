#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

log() {

        echo $(date -u)": "$1 >> $logfile

}

processGenome() {
   
    # default parameters
    local rel=2
    local k=9
    local E="1.0e-5"
    local path_to_db="/net/refdata/kegg"

    # parse arguments
    while [[ $# > 1 ]]
    do
    key=$1
        case $key in 
            -id | --genome_id)
                local genome_id=$2
                shift
                ;;
            -f | --assembly)
                local assembly=$2
                shift
                ;;
            -o | --output_dir)
                local output_dir=$2
                shift
                ;;
            -k | --kmer)
                local k=$2
                shift
                ;;
            -r | --reliability)
                local rel=$2
                shift
                ;;
            -o | --output_dir)
                local output_dir=$2
                shift
                ;;
            -E | --eval_threshold)
                local E=$2
                shift
                ;;
            -db | --path_to_db)
                local path_to_db=$2
                shift
                ;;
            -n | --n_cores)
                local n_cores=$2
                shift
                ;;
            *)
               # unknown option
               ;;
        esac
    shift
    done

    log "["$genome_id"] processing genome..."

    geneCalling -id $genome_id -f $assembly -o $output_dir
    prepareProteome -id $genome_id -ffn $output_dir/"$genome_id".ffn -faa $output_dir/"$genome_id".faa -o $output_dir
    annotateGenomeKEGG -id $genome_id -E $E -faa $output_dir/"$genome_id".faa -db $path_to_db -n $n_cores -o $output_dir

    log "["$genome_id"] done"

}

geneCalling() {

    # parse arguments
    while [[ $# > 1 ]]
    do
    key="$1"
        case $key in 
            -id | --genome_id)
                local genome_id="$2"
                shift
                ;;
            -f | --assembly)
                local assembly="$2"
                shift
                ;;
            -o | --output_dir)
                local output_dir="$2"
                shift
                ;;
            *)
               # unknown option
               ;;
        esac
    shift
    done

    log "["$genome_id"] gene calling..."

    # cleanup
    rm -f $output_dir/"genome_id".f*

    # use PRODIGAL for gene prediction
    prodigal -i $assembly -a $output_dir/"$genome_id".faa \
                          -d $output_dir/"$genome_id".ffn \
                          &>> $output

    # remove stop codons
    sed -i 's/\*//g' $output_dir/"$genome_id".faa

}

prepareProteome() {

    # parse arguments
    while [[ $# > 1 ]]
    do
    key="$1"
        case $key in 
            -id | --genome_id)
                local genome_id=$2
                shift
                ;;
            -ffn | --proteome_nt)
                local proteome_nt=$2
                shift
                ;;
            -faa | --proteome_aa)
                local proteome_aa=$2
                shift
                ;;
            -o | --output_dir)
                local output_dir=$2
                shift
                ;;
           *)
               # unknown option
               ;;
        esac
    shift
    done

    log "["$genome_id"] parsing proteome fasta files..."
  

    # cleanup
    rm -f $output_dir/"$genome_id"_compliant.ffn \
          $output_dir/"$genome_id"_compliant.faa

    # parse NT proteome
    awk -v org=$genome_id 'BEGIN {n=1}; />/ {print ">"org"|"org".peg."n; n++} !/>/ {print}' $proteome_nt \
        >> $output_dir/"$genome_id"_compliant.ffn
     mv $output_dir/"$genome_id"_compliant.ffn $output_dir/"$genome_id".ffn
    
    # parse AA proteome
    awk -v org=$genome_id 'BEGIN {n=1}; />/ {print ">"org"|"org".peg."n; n++} !/>/ {print}' $proteome_aa \
        >> $output_dir/"$genome_id"_compliant.faa
     mv $output_dir/"$genome_id"_compliant.faa $output_dir/"$genome_id".faa
    
}

annotateGenomeKEGG() {

    local path_to_db="/net/refdata/kegg"
    local n_cores=1
    local min_cov=90
    local E=0.001

    # parse arguments
    while [[ $# > 1 ]]
    do
    key="$1"
        case $key in 
            -id | --genome_id)
                local genome_id=$2
                shift
                ;;
            -faa | --proteome)
                local proteome=$2
                shift
                ;;
            -o | --output_dir)
                local output_dir=$2
                shift
                ;;
            -E | --eval_threshold)
                local E=$2
                shift
                ;;
            -c | --min_coverage)
                local min_cov=$2
                shift
                ;;
            -db | --path_to_db)
                local path_to_db=$2
                shift
                ;;
            -n | --n_cores)
                local n_cores=$2
                shift
                ;;
            *)
               # unknown option
               ;;
        esac
    shift
    done

    log "["$genome_id"] annotating genome using KEEG..."
    
    # cleanup
    rm -f $output_dir/"$genome_id".ko* $working_dir/"$genome_id"_ko*

    local kos_list=$(find $path_to_db/custom/hmms/ | grep -o 'K[0-9]*.hmm$' | sed 's/.hmm//g')
    
    for ko in $kos_list
    do
        
        hmmsearch --cpu $n_cores \
                  -E $E \
                  --domtblout $output_dir/"$genome_id"_ko_"$ko".txt \
                  $path_to_db/custom/hmms/"$ko".hmm $proteome \
                  &>> $output
       
        cat $output_dir/"$genome_id"_ko_"$ko".txt | \
            sed '/#/d;s/  */\t/g' | cut -f 1,3,4,6,7,16-19 | \
            awk -v min_c="$min_cov" 'BEGIN {
                                         IFS=FS=OFS="\t"
                                     }
                                     {
                                          hmm_c=($7-$6)*100/$4
                                          q_c=($9-$8)*100/$2
                                          if (q_c>min_c) {
                                               print $1,$3,$5,hmm_c,q_c
                                          }

                                     }' \
            >> $output_dir/"$genome_id"_ko_res.txt

        cat $output_dir/"$genome_id"_ko_res.txt | \
            awk -v E="$E" 'BEGIN {IFS=FS=OFS="\t"} {if ($3 < E) {print $0}}' \
            >> $output_dir/"$genome_id"_ko_res_filtered.txt
        
        mv $output_dir/"$genome_id"_ko_res_filtered.txt $output_dir/"$genome_id"_ko_res.txt

        rm -f $output_dir/"$genome_id"_ko_"$ko".txt
   
    done

    # parsing annotation results
   
    # cleanup 
    log "parsing results..."
    rm -f $output_dir/"$genome_id".ko_all \
        $output_dir/"$genome_id".ko \
        $output_dir/"$genome_id".ko_sorted \
        $output_dir/"$genome_id"_KEGG.txt

    # take best hit per peg ID

    pegs=$(grep ">" $proteome | sed 's/>//g')

    for i in $pegs
    do
        echo $i $(grep ''$i'[^0-9]' $output_dir/"$genome_id"_ko_res.txt | sort -g -k3 | cut -f 2 | head -n 1) \
            >> $output_dir/"$genome_id".ko
    done

    sed -i 's/ /\t/g' $output_dir/"$genome_id".ko

    # make final table (peg class pathway ko description)
   
    # sort peg ko table by ko term 
    sort -k 2,2  $output_dir/"$genome_id".ko >> $output_dir/"$genome_id".ko_sorted
    
    # join peg ko table with ko description
    join -t$'\t' $output_dir/"$genome_id".ko_sorted \
                 $path_to_db/custom/class-pathway-ko-description.txt \
                 -1 2 -2 4 -o 1.1,2.1,2.2,2.3,2.4,2.5 -a 1 | \
                 sed 's/[^\t]*peg.//g' | sort -t$'\t' -n -k1,1 | \
                 sed 's/^/'$genome_id'\|'$genome_id'.peg./g' \
                 >> $output_dir/"$genome_id"_KEGG.txt

    # cleanup
    rm -f $output_dir/"$genome_id".ko_all \
          $output_dir/"$genome_id".ko_sorted
          # $output_dir/"$genome_id"_ko_res.txt

}

