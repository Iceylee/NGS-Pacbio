samtools view -u file.sam | samtools sort -o file_sorted

stringtie $datadir/$i -p 10 -o stringtie_out/${foo}.gtf -G /data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf -l $foo

stringtie --merge -G $gtf -p 8 -o merged_EV.gtf mergelist_EV.txt




samtools view -u HBR_1.sam | samtools sort -o HBR_1.sorted.bam

rule stringtie:
	input:"outData/sorted_bam/{sample}_sorted.bam"
	output:"outData/stringtie/{sample}.gtf"
	params:gtf=config["ref"]["gtf"]
	shell:"""
             stringtie {input} -p 10 -o {output} -G {param.gtf} -l {wildcards.sample}


"""


rule mergelist:
        input:expand("outData/stringtie/{sample}.gtf",sample=SAMPLES)
        output:"outData/stringtie/mergelist_{condition}.txt"
        shell:"""

                        for file in `ls outData/stringtie/{wildcards.condition}_*.gtf`
                        do
                                echo $file >> {output}
                        done
"""

rule stringtie_merge:
        input:"outData/stringtie/mergelist_{condition}.txt"
        output:"outData/as/merge_{condition}.gtf"
        shell:"""
                     stringtie --merge -G config["ref"]["gtf"] -p 8 -o {output} {input}

"""

rule as:
		input:"outData/as/merge_{condition}.gtf"
		output:"outData/as/{condition}_event_count.txt"
		shell:"""
                 /data1/software/astalavista-4.0/bin/astalavista -t asta --threads 10 -i {input}
                 gunzip merge_{wildcards.condition}_sorted.gtf_astalavista.gtf.gz
                 awk -f scripts/as.awk merged_{wildcards.condition}_sorted.gtf_astalavista.gtf > {wildcards.condition}_as_count.txt
  

"""


/data1/software/astalavista-4.0/bin/astalavista -t asta --threads 10 -i outData/as/merge_UHR.gtf





