rule count_matrix:
    input:
        expand("outData/htseq/{sample}_CountNum.txt",sample=SAMPLES)
    output:
        "outData/counts/all.csv"
    params:
        units=units
    script:
        "../scripts/count-matrix.py"

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
    input:
        countTab="outData/counts/all.csv",
        colData=config["diffexp"]["colData"]
    output:
        rds="outData/deseq2/all.rds",
        countMatrix="outData/deseq2/norm_count.matrix.txt"
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule heatmap_cor_plot:
    input:
         colData=config["diffexp"]["colData"],
         countMatrix="outData/deseq2/norm_count.matrix.txt"
    output:
         heatmap_cor_plot=report("outData/deseq2/heatmap_cor.png")
    conda:
         "../envs/deseq2.yaml"
    log:
        "logs/deseq2/heatmap_cor_plot.log"
    script:
        "../scripts/heatmap_cor_plot.R"

rule heatmap_plot:
    input:
         colData=config["diffexp"]["colData"],
         rds="outData/deseq2/all.rds"
    output:
         heatmap_plot=report("outData/deseq2/heatmap.png")
    conda:
         "../envs/deseq2.yaml"
    log:
        "logs/deseq2/heatmap_plot.log"
    script:
        "../scripts/heatmap_plot.R"

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

rule deseq2:
    input:
        rds="outData/deseq2/all.rds"
    output:
        all_tab=report("outData/deseq2/{contrast}_all_genes_exprData.txt", "../report/diffexp.rst"),
        sig_tab=report("outData/deseq2/{contrast}_sig_genes_exprData.txt", "../report/diffexp.rst"),    
        plot=report("outData/deseq2/{contrast}_volcano_plot.png")
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/diffexp.R"
