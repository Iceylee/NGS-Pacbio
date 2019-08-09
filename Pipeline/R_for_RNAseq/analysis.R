#!/usr/bin/env Rscript

# scripts to complete differential expression and GO/KEGG enrichment analysis 
#
# originally by Yubing Li
# chashuguzi@163.com

# 1. Differential Expression Analysis===========================================
source("src/diff_express_fxns.R")

getNormalizeMat("data/raw/htseq_out.csv", "data/raw/colData.csv", "HBR_vs_UHR", "data/processed/norm_count_matrix.txt", "data/processed/dds.RData")

plotSmpClustHeatmap("data/processed/norm_count_matrix.txt", "data/raw/colData.csv", "output/heatmap_cor.pdf", "output/heatmap_cor.png", height = 6, width = 8, legend = TRUE,display_numbers = TRUE, number_format = "%.3f", number_color="black",cluster_cols = FALSE,cluster_rows = FALSE,annotation_col=anno,annotation_row=anno,annotation_colors = colors,fontsize_row = 12)

plotCountMatHeatmap("data/processed/dds.RData", 50, "output/count_heatmap.pdf", "output/count_heatmap.png", cluster_cols = FALSE,fontsize_row = 12)


getDiffTab("data/processed/dds.RData", "HBR_vs_UHR", "output/all_genes_exprData.txt", "output/sig_genes_exprData.txt", padj.threshold=1, pvalue.threshold=0.05, fold.threshold=1)

plotVolcanoGraph("data/processed/dds.RData", "HBR_vs_UHR", padj.threshold=0.05, fold.threshold=1, "output/volcano_plot.pdf", "output/volcano_plot.png")

# 2. GO & KEGG analysis=========================================================
Ah_id = getKEGGName("Homo sapiens")["Ah_id"]
KEGG_id = getKEGGName("Homo sapiens")["KEGG_id"]

mkdir(pset=1)
mkdir(pset=0.05)