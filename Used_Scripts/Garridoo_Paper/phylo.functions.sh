#!/bin/bash

# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

log() {

        echo $(date -u)": "$1 >> $logfile

}

getAmphoraGenes() {


    local genome_id=$1

    log "["$genome_id"] retrieving AMPHORA genes..."

    # cleanup
    rm -f $working_dir/"$genome_id"_amphora.txt

    local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

    for gene in $amphora_list
    do

        hmmsearch -E 0.001 \
                  -o $working_dir/"$genome_id"_"$gene".hmr \
                  $path_to_db/amphora/"$gene".hmm \
                  $data_dir/ORFs/"$genome_id".faa \
                  &> $output

        grep -o '....\.peg\.[0-9]*' $working_dir/"$genome_id"_"$gene".hmr | head -n 1 \
             >> $working_dir/"$genome_id"_"$gene".txt

        rm -f $working_dir/"$genome_id"_"$gene".hmr

        touch $working_dir/"$gene"_"$genome_id".ffn

        n_seqs=$(wc -l $working_dir/"$genome_id"_"$gene".txt | sed 's/ .*//g')

        if [[ "$n_seqs" -gt 0 ]]
        then

            peg=$(cat $working_dir/"$genome_id"_"$gene".txt)
            
            echo ">"$genome_id"" > $working_dir/"$gene"_"$genome_id".ffn
      	    awk "/$peg$/ {flag=1;next} />/{flag=0} flag {print}" $data_dir/ORFs/"$genome_id".ffn \
                >> $working_dir/"$gene"_"$genome_id".ffn

        fi

        # cleanup
        rm -f $working_dir/"$genome_id"_"$gene".txt

    done

}

alignAmphoraSeqs() {


    local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

    for gene in $amphora_list
    do
    
        log "["$gene"] alignning gene..."

        cat $working_dir/"$gene"_*.ffn >> $working_dir/"$gene".ffn
        rm -f $working_dir/"$gene"_*.ffn $working_dir/"$gene"*.msa

        clustalo --seqtype=DNA \
                 --threads=$(nproc) \
                 -i $working_dir/"$gene".ffn \
                 -o $working_dir/"$gene".msa \
                 --percent-id \
                 --distmat-out=$working_dir/"$gene"_distmat_raw.txt \
                 --full \
                 &>> $output
        
        tail -n +2 $working_dir/"$gene"_distmat_raw.txt >> $working_dir/"$gene"_distmat.txt
        rm -f $working_dir/"$gene"_distmat_raw.txt

    done


}

concatenateAmphoraSeqs() {


    rm -f $working_dir/amphora.msa $working_dir/amphora_list.txt

    local genome_list=$(cut -f 1 $data_dir/mapping.txt | grep -v "ID")
    local amphora_list=$(cat $path_to_db/amphora/amphora_list.txt)

    for genome_id in $genome_list 
    do

    echo ">"$genome_id"" >> $working_dir/amphora.msa

        for gene in $amphora_list
        do 
            
            local n=$(grep -c ">" $working_dir/"$gene".msa)
            local n_tot=$(echo $genome_list | awk '{print NF}')

            # use only genes found in every genome
            if [[ "$n" -eq $n_tot ]]
            then 

                awk "/$genome_id$/ {flag=1;next} />/{flag=0} flag {print}" $working_dir/"$gene".msa \
                    >> $working_dir/amphora.msa

                echo $gene >> $working_dir/amphora_list.txt

            fi

        done  

    done

}

buildBayesTree() {

    # default parameters
    local ngen=300000

    local ngen=$1

    log "building phylogeny..."

    rm -f $working_dir/amphora.nex

    clustalw -INFILE=$working_dir/amphora.msa \
             -OUTFILE=$working_dir/amphora.nex \
             -OUTPUT=NEXUS \
             -CONVERT \
             &>> $output

    cat mrbayes_block.nex | sed 's/NITER/'$ngen'/g' >> $working_dir/amphora.nex
    
    # run MrBayes and cleanup
    
    cd $working_dir
    mb $working_dir/amphora.nex &>> $output
    mv mrbayes_tmp.con.tre $working_dir/amphora_tree.tree
    rm -f $working_dir/mrbayes_tmp* $working_dir/amphora.nex
    cd $compgen_dir

}

buildMLTree() {

    # cleanup
    rm -f $working_dir/amphora.tree

    # generate ML ssp. tree
    FastTree -nt -gtr $working_dir/amphora.msa >> $working_dir/amphora.tree 2>> $output

}

