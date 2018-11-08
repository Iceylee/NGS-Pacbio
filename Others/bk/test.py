

hisat2-build ../refs/22.fa 22

rule hisat2:
	input:fwd="data/trimmed/{sample}_{unit}_clean_R1.fq.gz",
		  rev="data/trimmed/{sample}_{unit}_clean_R2.fq.gz"

	output:sam="data/hisat2/{sample}_{unit}.sam"
		   log="data/hisat2/{sample}_{unit}.log"

	params:index=config["ref"]["index"]
	message:"""--- Hisat2 Mapping.---"""
	shell:"""
 		 hisat2 -p 15 -x {params.index} -1 {input.fwd} -2 {input.rev} -S {output.sam} > {output.log} 2>&1

"""

rule sam_sort:
	input:"data/hisat2/{sample}_{unit}.sam"

	output:"data/sorted_bam/{sample}_{unit}_sorted.bam"
	shell:"""
		  samtools view -u {input} | samtools sort -@ 15 - > {output}
"""


htseq-count -f sam $samfile $gtf -q > $HtseqOut &

rule htseq:
	input:"data/hisat2/{sample}_{unit}.sam"

	output:"data/htseq/{sample}_{unit}_CountNum.txt"
    params:gtf=config["ref"]["gtf"]
	message:"""---htseq count---"""
	shell:"""
         htseq-count -f sam {input} {params.gtf} -q > {output} &


"""