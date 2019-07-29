rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        units=units
    script:
        "../scripts/count-matrix.py"

#####count-matrix.py###
import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in snakemake.input]

for t, (sample, unit) in zip(counts, snakemake.params.units.index):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum()
print(matrix)
matrix.to_csv(snakemake.output[0], sep="\t")
###########################

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


rule deseq2_init:
    input:
        countTab="data/counts/all.csv",
        colData=[config]["colData"]
    output:
        rds="data/deseq2/all.rds",
        countMatrix="data/deseq2/norm_count.matrix.txt"
    params:
        samples=config["samples"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"





rule heatmap_plot:
    input:
         rds="deseq2/all.rds"
    output:
         heatmap_plot=report("results/diffexp/heatmap.pdf")
    conda:
         "../envs/deseq2.yaml"
    log:
        "logs/deseq2/heatmap_plot.log"
    script:
        "../scripts/heatmap_plot.R"

rule heatmap_cor_plot:
    input:
         colData=[config]["colData"]
         countMatrix="data/deseq2/norm_count.matrix.txt"
    output:
         heatmap_cor_plot=report("results/diffexp/heatmap_cor.pdf")
    conda:
         "../envs/deseq2.yaml"
    log:
        "logs/deseq2/heatmap_cor_plot.log"
    script:
        "../scripts/heatmap_cor_plot.R"






def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

rule deseq2:
    input:
        "data/deseq2/all.rds"
    output:
        all_tab=report("data/deseq2/{contrast}_all_genes_exprData.txt", "../report/diffexp.rst"),
        sig_tab=report("data/deseq2/{contrast}_sig_genes_exprData.txt", "../report/diffexp.rst")
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/diffexp.R"


