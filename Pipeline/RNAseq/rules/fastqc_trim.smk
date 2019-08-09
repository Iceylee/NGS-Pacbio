def get_fastq(wildcards):
    return units.loc[wildcards.sample, ["fq1", "fq2"]].dropna()

rule trimming:
  input: get_fastq
  output: fwd="outData/trimmed/{sample}_clean_R1.fq.gz",
          rev="outData/trimmed/{sample}_clean_R2.fq.gz",
	  up_fwd="outData/trimmed/unpaired_{sample}_R1.fq.gz",
	  up_rev="outData/trimmed/unpaired_{sample}_R2.fq.gz"
  threads: config["threads"]
  message: """--- Trimming."""
  shell: """
         java -jar /data1/software/Trimmomatic-0.36/trimmomatic-0.36.jar  PE -threads {threads}  -phred33  {input} {output.fwd} {output.up_fwd} {output.rev} {output.up_rev} ILLUMINACLIP:/data1/software/Trimmomatic-0.36/adapters/TruSeq3-PE_all.fa:2:30:10:8:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:18 MINLEN:36 >> trimmomatic.log
         #rm {output.up_fwd} {output.up_rev}
  """
rule fastqc:
        input: get_fastq 
        output: fwd="outData/fastqc/{sample}_R1_fastqc.zip",
                rev="outData/fastqc/{sample}_R2_fastqc.zip"
        message: """--- Quality check of raw outData with Fastqc."""
        shell: """
              fastqc --outdir outData/fastqc/  --extract  -f fastq {input}
     """
rule fastqc_clean:
        input:fwd="outData/trimmed/{sample}_clean_R1.fq.gz",
              rev="outData/trimmed/{sample}_clean_R2.fq.gz"
        output: fwd="outData/fastqc_clean/{sample}_clean_R1_fastqc.zip",
                rev="outData/fastqc_clean/{sample}_clean_R2_fastqc.zip"
        message: """--- Quality check of clean outData with Fastqc."""
        shell: """
              fastqc --outdir outData/fastqc_clean/  --extract  -f fastq {input}
     """
