#!/usr/bin/env Rscript

arg <- commandArgs(T)

# test argument
if (length(args) != 2) {
	print ("2 arguments must be supplied: ")
	print ("outputPath otu_table_mc2_w_tax_no_pynast_failures.txt L6.txt")
	stop(call.=FALSE)
}

library(pheatmap)

outputPath <- argv[1]
otu_file <- argv[2]
otu_L6_file <- argv[3]

correlation_tif <- paste(outputPath, "8.Other_Analysis/correlation.tif", sep="")

otu_table <- read.delim(otu_file,header=TRUE,row.names=1)
otu_table_colsorted <- otu_table[,order(names(otu_table))]
cortable <- cor(otu_table_colsorted)

# 1 heatmap
tiff(correlation_tif,height=6,width=8,units='in',res=300)
pheatmap(cortable,cluster_rows = FALSE,cluster_cols = FALSE, border_color = "grey95",scale = "none",display_numbers = TRUE, number_format = "%.2f", number_color="black")
dev.off()


taxa_table <- read.delim(otu_L6_file,header=TRUE,row.names=1)
taxa_table[taxa_table==0] <- 0.0000000001
taxa_table0.01 <- taxa_table[rowSums(taxa_table>0.01)>=1,]

OTU_Heatmap_tif <- paste(outputPath, "8.Other_Analysis/OTU_heatmap.tif", sep="")

# 2 Heatmap
#high_taxa_table <- taxa_table[head(order(rowSums(taxa_table),decreasing=TRUE),30),]
tiff(OTU_Heatmap_tif,height=10,width=15,units='in',res=300) 
pheatmap(mat=log10(taxa_table0.01_colsorted),scale="row",fontsize_row=9)
dev.off()

#write.taxa_table(taxa_table0.01,file="key_OTU_relative_abundance.txt", quote=FALSE, sep="\t",col.names=NA)


