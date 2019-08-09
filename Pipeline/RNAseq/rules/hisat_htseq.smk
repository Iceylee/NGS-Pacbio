rule hisat2:
	input:fwd="outData/trimmed/{sample}_clean_R1.fq.gz",
		  rev="outData/trimmed/{sample}_clean_R2.fq.gz"
	output:sam="outData/hisat2/{sample}.sam",
	       log="outData/hisat2/{sample}.log"
	params:index=config["ref"]["index"]
        threads:config["threads"]
	message:"""--- Hisat2 Mapping.---"""
	shell:"""
 		 hisat2 -p {threads} -x {params.index} -1 {input.fwd} -2 {input.rev} -S {output.sam} > {output.log} 2>&1

"""

rule sam_sort:
	input:"outData/hisat2/{sample}.sam"
	output:"outData/sorted_bam/{sample}_sorted.bam"
        threads:config["threads"]
	shell:"""
		  samtools view -u {input} | samtools sort -@ {threads} - > {output}
"""

rule htseq:
	input:"outData/hisat2/{sample}.sam"
	output:"outData/htseq/{sample}_CountNum.txt"
        params:gtf=config["ref"]["gtf"]
	message:"""---htseq count---"""
	shell:"""
         htseq-count -f sam {input} {params.gtf} -q > {output} 

"""
