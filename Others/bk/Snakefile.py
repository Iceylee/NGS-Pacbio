from snakemake.utils import R
import sys

"""
Author: D. Puthier
Affiliation: AMU
Aim: A simple Snakemake workflow to process paired-end stranded RNA-Seq.
Date: Mon Nov  2 14:03:11 CET 2015
Run: snakemake   -s Snakefile   
Latest modification: 
  - todo
"""

##-----------------------------------------------##
## A set of functions
##-----------------------------------------------##

def message(mes):
  sys.stderr.write("|--- " + mes + "\n")

##-----------------------------------------------##
## Working directory                             ##
## Adapt to your needs                           ##
##-----------------------------------------------##

BASE_DIR="/Users/puthier/git/abd/practical/snakemake/projects/"
WDIR = BASE_DIR + "t-cell_analysis"
workdir: WDIR
message("The current working directory is " + WDIR)

##--------------------------------------------------------------------------------------##
## Variables declaration                          
## Declaring some variables used by topHat and other tools... 
## (GTF file, INDEX, chromosome length)
##--------------------------------------------------------------------------------------##
# Adapt the path to your needs
INDEX = BASE_DIR + "indexes/bowtie2/mm10/chr19"
GTF   = BASE_DIR + "annotations/gtf/GRCm38.83.chr19.gtf"
CHR   = BASE_DIR + "annotations/txt/chromInfo_mm10.txt"
FASTA = BASE_DIR + "indexes/bowtie2/mm10/chr19.fa"

##--------------------------------------------------------------------------------------##
## The list of samples to be processed
##--------------------------------------------------------------------------------------##

SAMPLES, = glob_wildcards("samples/raw/{smp}_R1.fq.gz")

NB_SAMPLES = len(SAMPLES)


for smp in SAMPLES:
  message("Sample " + smp + " will be processed")

rule final: 
  input:  expand("samples/fastqc/{smp}/{smp}_R1_t.fq_fastqc.zip", smp=SAMPLES),\
          expand("samples/bam/{smp}.bam.bai", smp=SAMPLES), \
          expand("samples/bwig/{smp}_min.bw", smp=SAMPLES), \
          "results/diagnostic_plot/diagnostic.pdf"

rule trimming:
  input:  fwd="samples/raw/{smp}_R1.fq.gz", rev="samples/raw/{smp}_R2.fq.gz"
  output: fwd="samples/trimmed/{smp}_R1_t.fq", 
          rev="samples/trimmed/{smp}_R2_t.fq", 
          single="samples/trimmed/{smp}_R1_singletons.fq"
  message: """--- Trimming reads."""
  shell: """                                                                                                                                                                                                                                                                                                                                                            
    sickle pe -f {input.fwd} -r {input.rev}  -l 40 -q 20 -t sanger  -o {output.fwd} -p {output.rev} -s {output.single} &> {input.fwd}.log                                                                               
  """

rule fastqc:
        input:  fwd="samples/trimmed/{smp}_R1_t.fq", 
                rev="samples/trimmed/{smp}_R2_t.fq"
        output: fwd="samples/fastqc/{smp}/{smp}_R1_t.fq_fastqc.zip", rev="samples/fastqc/{smp}/{smp}_R2_t.fq_fastqc.zip"
        message: """--- Quality check of raw data with Fastqc."""
        shell: """                                                                                                                                                                                        
        /opt/FastQC/fastqc --outdir  samples/fastqc/{wildcards.smp} --extract  -f fastq {input.fwd} {input.rev}    &> {output}.log                                                                                                                                                                                                                                                                                
     """

rule tophat:
        input:fwd="samples/trimmed/{smp}_R1_t.fq", 
              rev="samples/trimmed/{smp}_R2_t.fq"
        params: gtf=GTF, index=INDEX
        output: "samples/bam/{smp}.bam"
        message: "--- Mapping reads with tophat."
        shell: """
                mkdir -p samples/tophat/{wildcards.smp}
                tophat2                                         \
                        -o samples/tophat/{wildcards.smp}    \
                        -g 1                                    \
                        --library-type fr-firststrand           \
                        -G {params.gtf}                          \
                        -x 1                                    \
                        -p 5                                    \
                        {params.index}                             \
                        {input.fwd} {input.rev} &> samples/tophat/{wildcards.smp}/run_tophat.log
                cd samples/tophat/{wildcards.smp}
                mv accepted_hits.bam ../../bam/{wildcards.smp}.bam
        """

rule select_reads_by_strand:
        input: "samples/bam/{smp}.bam"
        output: min="samples/bam/{smp}_min.bam", plus="samples/bam/{smp}_plus.bam"
        shell: """
              samtools  view -f99 -hb {input} > {output.min}_99.bam
              samtools  view -f147 -hb {input} > {output.min}_147.bam
              samtools  view -f83 -hb {input} > {output.plus}_83.bam
              samtools  view -f163 -hb {input} > {output.plus}_163.bam
      
              samtools merge -f -h {output.min}_99.bam  {output.min} {output.min}_99.bam {output.min}_147.bam
              samtools merge -f -h {output.plus}_83.bam {output.plus} {output.plus}_83.bam {output.plus}_163.bam
      
              rm -f {output.min}_99.bam  {output.min}_147.bam {output.plus}_83.bam {output.plus}_163.bam
        """

rule index_bam:
        input: full="samples/bam/{smp}.bam", min="samples/bam/{smp}_min.bam", plus="samples/bam/{smp}_plus.bam"
        output: full="samples/bam/{smp}.bam.bai"
        message: "--- Indexing BAM files."
        shell: """
              samtools index {input.full}
              samtools index {input.min}
              samtools index {input.plus}
        """

rule create_bigwig:
        input: min="samples/bam/{smp}_min.bam", plus="samples/bam/{smp}_plus.bam"
        output: min="samples/bwig/{smp}_min.bw", plus="samples/bwig/{smp}_plus.bw"
        params: chr=CHR
        message: "--- Creating bigwig from BAM files."
        shell: """
          bam2wig.py  -i {input.plus} -s {params.chr} -o samples/bwig/{wildcards.smp}_plus &> {output.plus}.log
          wigToBigWig -clip samples/bwig/{wildcards.smp}_plus.wig {params.chr} {output.plus} 2>&1 >> {output.plus}.log
          rm -f samples/bwig/{wildcards.smp}_plus.wig
  
          bam2wig.py  -i {input.min} -s {params.chr} -o samples/bwig/{wildcards.smp}_min  &> {output.min}.log
          wigToBigWig -clip samples/bwig/{wildcards.smp}_min.wig {params.chr} {output.min} 2>&1 >> {output.min}.log
          rm -f samples/bwig/{wildcards.smp}_min.wig
      """

rule cufflinks:
        input: bam="samples/bam/{smp}.bam"
        output: gtf="samples/cufflinks/{smp}/transcripts.gtf"
        params:  gtf=GTF
        message: "--- Searching novel transcript with cufflinks."
        shell: """
          cufflinks -g {params.gtf} -p 5 --library-type fr-firststrand  -o samples/cufflinks/{wildcards.smp}  {input.bam} &> {output}.log
        """

rule cuffmerge:
        input: expand("samples/cufflinks/{smp}/transcripts.gtf", smp=SAMPLES)
        output: "samples/cuffmerge/merged.gtf"
        params: gtf=GTF, fa=FASTA
        message: "--- Comparing transcript to the reference."
        shell:  """
        ls -1 samples/cufflinks/*/transcripts.gtf > samples/cuffmerge/assembly.txt
        cuffmerge  -o  samples/cuffmerge -g {params.gtf} --keep-tmp -s {params.fa} -p 5 samples/cuffmerge/assembly.txt &> {output}.log
        """

rule select_novel_transcripts:
        input: "samples/cuffmerge/merged.gtf"
        output: "samples/cuffmerge/novel_transcript.gtf"
        params: gtf=GTF, fa=FASTA
        message: "--- Selecting novel_transcripts."
        shell:  """
        awk '/class_code "[uxosj]"/' {input} > {output}
        """

rule add_gene_name_to_unknown:
        input: "samples/cuffmerge/novel_transcript.gtf"
        output: "samples/cuffmerge/novel_transcript_gn.gtf"
        params: gtf=GTF, fa=FASTA
        message: "--- Adding gene name to novel transcript."
        run:
          import re 
          fh_in = open(input[0], "r")
          fh_out = open(output[0], "w")
          for line in fh_in:
            line = line.rstrip("\n")
            if not re.search("gene_name", line):
              gene_id = re.match('.*gene_id "(.*?)"', line).group(1)
              fh_out.write(line + ' gene_name "' + gene_id + '";\n')

rule merge_novel_and_known:
        input: novel="samples/cuffmerge/novel_transcript_gn.gtf", known=GTF
        output: "samples/new_annotation/all_transcripts.gtf"
        params: gtf=GTF, fa=FASTA
        message: "--- Merging known and novel transcripts."
        shell:  """
        cat {input.novel} {input.known} > {output}
        """ 



rule quantification_with_featureCounts:
        input: novel="samples/new_annotation/all_transcripts.gtf", bam=expand("samples/bam/{smp}.bam", smp=SAMPLES)
        output: "results/counts/gene_counts.txt",  "results/counts/gene_counts_mini.txt"
        shell: """
        featureCounts -p -s 2 -T 15 -t exon -g gene_id -a {input.novel} -o {output[0]} {input.bam} &> {output[0]}.log
        cut -f 1,7- {output[0]}| awk 'NR > 1' | awk '{{gsub("samples/bam/","",$0); print}}'  > {output[1]}
        """

rule diagnostic_plot:
        input: "results/counts/gene_counts_mini.txt"
        output: "results/diagnostic_plot/diagnostic.pdf"
        run: R("""
            dir.create("results/diagnostic_plot")
            data <- read.table("{input}", 
                                sep="\t", 
                                header=T, 
                                row.names=1)
            data <- data[rowSums(data) > 0, ]
            data <- log2(data + 1)
            pdf("{output}")
            dev.null <- apply(data, 2, hist, border="white", col="blue")
            boxplot(data, color="blue", pch=16)
            pairs(data, pch=".", col="blue")
            dev.off()
            cat("etc...")
      
        """)


