rule stringtie:
	input:"outData/sorted_bam/{sample}_sorted.bam"
	output:"outData/stringtie/{sample}.gtf"
        params:gtf=config["ref"]["gtf"]
        threads:config["threads"]
	shell:"""
             stringtie {input} -p {threads} -o {output} -G {params.gtf} -l {wildcards.sample}


"""

rule mergelist:
    input:config["diffexp"]["colData"]
    output:"outData/stringtie/mergelist_{condition}.txt"
    params:con="{condition}"
    run:
        colData = pd.read_csv(input[0], dtype=str).set_index("sample", drop=False)
        SAMPLES = colData[colData.condition==params[0]].ix[:,'sample']
        output_file = open(output[0],"w")
        for i in SAMPLES:
            output_file.write("outData/stringtie/"+ i +".gtf"+"\n")
        output_file.close()

rule stringtie_merge:
        input:"outData/stringtie/mergelist_{condition}.txt"
        output:"outData/as/merge_{condition}.gtf"
        params:gtf=config["ref"]["gtf"]
        threads:config["threads"]
        shell:"""
                     stringtie --merge -G {params.gtf} -p {threads} -o {output} {input}

"""

rule as:
		input:"outData/as/merge_{condition}.gtf"
		output:"outData/as/{condition}_event_count.txt"
                threads:config["threads"]
		shell:"""
                 /data1/software/astalavista-4.0/bin/astalavista -t asta --threads {threads} -i {input}
                 gunzip outData/as/merge_{wildcards.condition}_sorted.gtf_astalavista.gtf.gz
                 awk -f scripts/as.awk outData/as/merge_{wildcards.condition}_sorted.gtf_astalavista.gtf > {output}
"""
