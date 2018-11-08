rule GATK_SelectVariants:
        input:colData=config["diffexp"]["colData"]
        output:"outData/snp/{condition}_merge.vcf"
        params:con="{condition}",fa=config["ref"]["fasta"]
        run:
            import os
            refName=params[1][:-3]
            cmd1 = "java -jar /data1/software/picard/picard.jar CreateSequenceDictionary R=%s O=%s.dict" % (params[1],refName)
            if !os.path.isfile("%s.dict" % (refName)):
                shell(cmd1)
            colData = pd.read_csv(input[0], dtype=str).set_index("sample", drop=False)
            SAMPLES = colData[colData.condition==params[0]].ix[:,'sample']
            temp=SAMPLES[0]
            for i in SAMPLES[1:]:
                temp=SAMPLES[0]
                cmd2 = "java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R %s --variant outData/snp/%s.samtools.vcf --concordance outData/snp/%s.samtools.vcf -o %s" % (params[1],temp,i,output[0])
                shell(cmd2)
                temp=output[0]



rule snpEff:
        input:"outData/snp/{condition}_merge.vcf"
        output:html="outData/snp/Annotation/{condition}.SnpEff.html",
               csv="outData/snp/Annotation/{condition}.SnpEff.csv"
        shell:"""
                    nohup java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s {output.html} -c snpEff.config  -v -ud 500 Sitotroga {input} > {output.csv} &
        """