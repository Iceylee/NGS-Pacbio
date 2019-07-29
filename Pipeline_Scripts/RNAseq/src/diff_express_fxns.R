# Load the library.
suppressMessages(library(DESeq2))
suppressMessages(library(dplyr))
suppressMessages(library(pheatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(genefilter))
suppressMessages(library(ggplot2))


#' Get deseq2 norm matrix.
#'
#' @param groups groupA_vs_groupB, group name should be consistent with 
#' coldata_file
#' @return
#' @example
#' getNormalizeMat("CountMatrix4DESeq.csv", "col_data.csv", "A_vs_B", 
#' "norm-count-matrix.txt", "dds.RData")
#' @author Yubing Li


getNormalizeMat <- function(count_file, coldata_file, groups, 
                            output_count_file, output_dds_file) {
  count_data <- read.csv(count_file, header = T, row.names = 1)
  col_data <- read.csv(coldata_file, header = T, row.names = 1)
  count_data <- count_data[, rownames(col_data)]
  col_data$condition <- as.factor(col_data$condition)

  exp_group <- strsplit(groups, "_vs_")[[1]][1]
  base_group <- strsplit(groups, "_vs_")[[1]][2]

  if (all(rownames(col_data) %in% colnames(count_data)) == F) {
    print("rownames(col_data) != colnames(count_data)")
    print("Please check.")
    stop(call. = FALSE)
  }

  # Create DESEq2 dataset.
  # ~在R里面用于构建公式对象，~左边为因变量，右边为自变量。
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, 
                                design = ~condition)

  # pre filter
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  # Set the reference to be compared
  dds$condition <- relevel(dds$condition, base_group)

  # Run deseq
  dds <- DESeq(dds)
  # Get normalized counts and write this to a file
  nc <- counts(dds, normalized = TRUE)
  # Turn it into a dataframe to have proper column names.
  dt <- data.frame("id" = rownames(nc), nc)
  colnames(dt) <- c("id", colnames(nc))
  # Save the normalized data matrix.
  write.table(dt, file = output_count_file, sep = "\t", row.name = FALSE, 
              col.names = TRUE, quote = FALSE)

  save(dds, file = output_dds_file)
}

#' Plot the sample clustering heatmap.
#'
#' @param
#' @return
#' @example
#' plotSmpClustHeatmap("norm-count-matrix.txt", "col_data.csv", 
#' "heatmap_cor.pdf", "heatmap_cor.png")
#' @author Yubing Li


plotSmpClustHeatmap <- function(norm_count_file, coldata_file, output_pdf, 
                                output_png, height = 6, width = 8, ...) {
  df <- read.table(file = norm_count_file, sep = "	", header = T, 
                   stringsAsFactors = F, row.names = 1)
  col_data <- read.csv(coldata_file, header = T, row.names = 1)
  colnames(df) <- rownames(col_data)
  cor.df <- cor(df[, unlist(lapply(df, is.numeric))])

  anno <- data.frame(Condition = col_data$condition)
  anno$Condition <- as.factor(anno$Condition)
  rownames(anno) <- rownames(col_data)

  # color sets
  colors <- c(brewer.pal(name = "Set1", n = 8), brewer.pal(name = "Paired", 
                                                           n = 12))
  colors <- colors[1:length(unique(anno$Condition))]
  names(colors) <- unique(anno[, 1])
  colors <- list(Condition = colors)

  # pdf plot
  pdf(file = output_pdf, height = height, width = width, onefile = FALSE)
  pheatmap(cor.df, annotation_col = anno, annotation_row = anno, 
           annotation_colors = colors, ...)
  dev.off()

  # png plot
  png(file = output_png, height = height, width = width, units = "in", 
      res = 300)
  pheatmap(cor.df, annotation_col = anno, annotation_row = anno, 
           annotation_colors = colors, ...)
  dev.off()
}


#' Plot the count matrix heatmap.
#'
#' @param
#' @return
#' @example
#' plotCountMatHeatmap("data/processed/dds.RData", 50, "output/count_heatmap.pdf", 
#' "output/count_heatmap.png", cluster_cols = FALSE,fontsize_row = 12)
#' @author Yubing Li

plotCountMatHeatmap <- function(dds_file, top_x, output_pdf, output_png,
                                height = 6, width = 8, ...) {
  load(dds_file)
  rld <- rlog(dds, blind = FALSE)

  # 前top_x样本差异较大的，做聚类热图。不同样本中的表达情况
  topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), top_x)
  mat <- assay(rld)[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  anno <- data.frame(Condition = col_data$condition)
  anno$Condition <- as.factor(anno$Condition)
  rownames(anno) <- rownames(col_data)

  # color sets
  # colors <- c("navy","darkgreen","yellow","red","blue","orange","grey",
  # "Blue Violet","Brown","Cyan")
  colors <- c(brewer.pal(name = "Set1", n = 8), brewer.pal(name = "Paired", 
                                                           n = 12))
  colors <- colors[1:length(unique(anno$Condition))]
  names(colors) <- unique(anno[, 1])
  colors <- list(Condition = colors)


  # pdf plot
  pdf(
    file = output_pdf, height = height, width = width,
    onefile = FALSE
  )
  pheatmap(mat, annotation_col = anno, annotation_colors = colors, ...)
  dev.off()

  # png plot
  png(
    file = output_png, height = height,
    width = width, units = "in", res = 300
  )
  pheatmap(mat, annotation_col = anno, annotation_colors = colors, ...)
  dev.off()
}


#' Differential analysis using DeSeq2.Get all genes' expression data and 
#' significant expression data.(Note: exlude padj NA genes)
#'
#' @param groups strings as "groupA_vs_groupB", group name should be consistent 
#' with coldata_file
#' @return
#' @example
#' getDiffTab("data/processed/dds.RData", "HBR_vs_UHR", 
#' "output/all_genes_exprData.txt", "output/sig_genes_exprData.txt", 
#' padj.threshold=1, pvalue.threshold=0.05, fold.threshold=1)
#' @author Yubing Li

getDiffTab <- function(dds_file, groups, output_diff_tab, output_sig_tab, 
                       padj.threshold = 0.05, pvalue.threshold = 0.05, 
                       fold.threshold = 1) {
  load(dds_file)

  exp_group <- strsplit(groups, "_vs_")[[1]][1]
  base_group <- strsplit(groups, "_vs_")[[1]][2]

  # 差异分析结果 multiple use contrast
  res <- results(dds, contrast = c("condition", exp_group, base_group))
  # Sort the results data frame by the padj and foldChange columns.
  sorted <- res[with(res, order(pvalue, -log2FoldChange)), ]
  # Turn it into a dataframe to have proper column names.
  sorted.df <- data.frame("id" = rownames(sorted), sorted)
  # Write the table out.
  write.table(sorted.df, file = output_diff_tab, row.names = FALSE, 
              sep = "\t", quote = FALSE)

  # significantly different genes
  res <- subset(res, padj <= padj.threshold)
  res <- subset(res, pvalue <= pvalue.threshold)
  res <- subset(res, abs(log2FoldChange) >= fold.threshold)
  res <- res[with(res, order(-log2FoldChange)), ]
  sig <- data.frame("id" = rownames(res), res)
  write.table(sig, file = output_sig_tab, sep = "\t", row.name = FALSE, 
              col.names = TRUE, quote = FALSE)
}

#' Plot volcano pdf and png for differential expression genes.
#'
#' @param groups strings as "groupA_vs_groupB", group name should be consistent 
#' with coldata_file
#' @return
#' @example
#' plotVolcanoGraph("data/processed/dds.RData", "HBR_vs_UHR", padj.threshold=0.05,
#'  fold.threshold=1, "output/volcano_plot.pdf", "output/volcano_plot.png")
#' @author Yubing Li
plotVolcanoGraph <- function(dds_file, groups, padj.threshold = 0.05, 
                             fold.threshold = 1, output_pdf, output_png) {
  exp_group <- strsplit(groups, "_vs_")[[1]][1]
  base_group <- strsplit(groups, "_vs_")[[1]][2]

  res <- results(dds, contrast = c("condition", exp_group, base_group))

  df <- data.frame("id" = rownames(res), res)
  sig_df <- filter(df, !is.na(pvalue))

  # scale_color_manual
  color <- c(red = "red", gray = "gray", blue = "blue")
  # add color column ,and condition
  sig_df$color <- ifelse(
    sig_df$padj < padj.threshold & abs(sig_df$log2FoldChange) >= fold.threshold, 
    ifelse(sig_df$log2FoldChange > fold.threshold, "red", "blue"), "gray")

  p <- ggplot(sig_df, aes(x = log2FoldChange, y = -log10(padj), col = color)) +
    geom_point() +
    scale_color_manual(values = color) +
    labs(x = "log2 (fold change)", y = "-log10 (padj)") +
    geom_hline(yintercept = -log10(0.05), lty = 4, col = "grey", lwd = 0.6) +
    geom_vline(xintercept = c(-1, 1), lty = 4, col = "grey", lwd = 0.6) +
    theme(
      legend.position = "none",
      panel.grid = element_blank()
    ) +
    theme(axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14)) +
    labs(title = paste(groups, " volcano plot", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5))

  ggsave(file = output_pdf, p, width = 6, height = 6, units = "in")

  ggsave(file = output_png, p, width = 6, height = 6, units = "in")
}