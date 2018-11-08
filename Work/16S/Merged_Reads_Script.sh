for file in $(<list.txt)
do
    join_paired_ends.py -f ${file}F_sub.fastq -r ${file}R_sub.fastq -o ${file}/
    mv ${file}/fastqjoin.join.fastq Merged_Reads/${file}_merged.fastq
    convert_fastaqual_fastq.py -c fastq_to_fastaqual -f Merged_Reads/${file}_merged.fastq -o Merged_Reads/${file}
    mv Merged_Reads/${file}/${file}_merged.fna Merged_Reads/${file}_merged.fasta
    rm -r Merged_Reads/${file}
    rm Merged_Reads/${file}_merged.fastq
    rm -r ${file}/
done
