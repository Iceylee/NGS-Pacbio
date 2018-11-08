
# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list = ls())

# load libraries

library(utils, quietly=T, warn.conflicts=F)
library(Biostrings, quietly=T, warn.conflicts=F)
library(grid, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)

options(warn=-1)

# plotting stuff

source("plotting.functions.R")

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line=element_line(color="black"),
              axis.ticks=element_line(color="black"),
              axis.text=element_text(colour="black", size=10),
              legend.position="top",
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(family="sans"))

# input data files

# load paths to project directories
source("paths.R")

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "wgs_taxonomy.txt", sep="")
amphora_list.file <- paste(amphora.dir, "amphora_list.txt", sep="")
functional.file <- paste(data.dir, "functional_profiles.txt", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")
tax <- read.table(taxonomy.file, sep="\t", header=T)
func <- read.table(functional.file, sep="\t", header=T) 
amphora_list <- read.table(amphora_list.file, sep="\t", header=F)[, 1] 

colnames(func) <- gsub("^X", "", colnames(func))

### functional MDS

message("calculating functional PCoA...")

d <- 1 - cor(func)
d[d==0] <- .1
diag(d) <- 0

k <- 2

pcoa <- cmdscale(d, k=k, eig=T)
mds <- pcoa$points
eig <- pcoa$eig
                                                                                
mds <- as.data.frame(mds)
colnames(mds) <- c("x", "y")
                                                                                
mds$compartment <- mapping$compartment[match(rownames(mds), mapping$ID)] 

tax$genome <- mapping$ID[match(tax$isolate_ID, mapping$Strain)]
mds$taxonomy <- tax$Phylum[match(rownames(mds), tax$genome)] 

### clustering isolates by family

tree  <- read.tree(paste(amphora.dir, "amphora_rooted_names.tree", sep=""))

groups_names <- rev(as.character(tax$Family[match(tree$tip.label, tax$isolate_ID)]))
groups_names <- match(groups_names, unique(groups_names))
names(groups_names) <- rev(tree$tip.label)
groups_ids <- groups_names
names(groups_ids) <- mapping$ID[match(names(groups_names), mapping$Strain)]

write.table(data.frame(ID=names(groups_names), cluster=groups_names),
            file=paste(amphora.dir, "/clusters.txt", sep=""),
            col.names=T, row.names=F, sep="\t", quote=F)

# read amphora MSA
amphora.msa <- readDNAStringSet(paste(amphora.dir, "amphora.msa", sep=""))

### functions needed to generate the plots

genTree <- function(k) {

    # subset sequences from the amphora MSA

    seqs <- amphora.msa[names(amphora.msa) %in% genomes]
    writeXStringSet(seqs, paste(amphora.dir, "amphora_tmp.msa", sep=""))

    # build sub-tree

    system(paste("rm -f ", amphora.dir, "amphora_tmp.tree", sep=""))

    system(paste("FastTree -nt -gtr ", amphora.dir, "amphora_tmp.msa >> ",
                 amphora.dir, "amphora_tmp.tree 2>> /dev/null", sep=""))

    tree  <- read.tree(paste(amphora.dir, "amphora_tmp.tree", sep=""))

    # calculate pair-wise distances between all leaves of the tree
    
    d <- cophenetic(tree)
    
    compartment <- mapping$compartment[match(colnames(d), mapping$ID)]
    
    idx.leaf <- compartment=="leaf"
    d.leaf <- dist(d[idx.leaf, idx.leaf])
    idx.root <- compartment=="root" | compartment=="soil"
    d.root <- dist(d[idx.root, idx.root])
    
    d <- dist(d)

    df.all <- data.frame(distances=as.vector(d), compartment="all")
    df.leaf <- if (sum(idx.leaf) > 1) data.frame(distances=as.vector(d.leaf), compartment="leaf") else NULL
    df.root <- if (sum(idx.root) > 1) data.frame(distances=as.vector(d.root), compartment="root") else NULL
    df <- rbind(df.all, df.leaf, df.root)

    # plot unrooted tree with color labels

    compartment <- mapping$compartment[match(tree$tip.label, mapping$ID)]

    tips.col <- rep(NA, length(tree$tip.label))
    tips.col[compartment=="leaf"] <- c_dark_green
    tips.col[compartment=="root"] <- c_red
    tips.col[compartment=="soil"] <- c_red

    tips.shape <- rep(NA, length(tree$tip.label))
    tips.shape[compartment=="leaf"] <- 24
    tips.shape[compartment=="root"] <- 21
    tips.shape[compartment=="soil"] <- 22

    tips.cex <- rep(1, length(tree$tip.label))
 
    tree$tip.label <- mapping$Strain[match(tree$tip.label, mapping$ID)]

    plot_btree(tree, output.file=paste(figures.dir, "/amphora_tmp_tree.pdf", sep=""),
               cex.labels=.7, label.offset=.03,
 	           type="unrooted", use.edge.length=T, edge.width=1,
               align.labels=T, color.labels="transparent", color.alg.lines="dark grey",
               tips.col=tips.col, tips.shape=tips.shape, tips.cex=tips.cex,
               margins=c(.1, 0, .2, 8), height=6.5, width=11.7)

    return(df)

}

genDensities <- function(k, df, points) {

    # plot boxplots

    colors1 <- data.frame(compartment=c("all", "leaf", "root"),
                          cols=c("black", c_dark_green, c_red))
    c1 <- as.character(colors1$cols[match(levels(df$compartment), colors1$compartment)])
    
    p1 <- ggplot(df, aes(x=compartment, y=distances, color=compartment, shape=compartment)) +
                 geom_boxplot(alpha=1, outlier.size=.5, size=.7) +
                 geom_jitter(position=position_jitter(0.35), size=.7, alpha=.25) +
                 scale_shape_manual(values=c(1, 17, 16)) +
                 scale_colour_manual(values=c1) +
                 labs(x="", y="mutations per site", title="pair-wise tree distances") +
                 coord_flip() +
                 main_theme +
                 theme(legend.title=element_blank())

    # plot MDS sub-plot

    points$compartment <- as.character(points$compartment)
    points$compartment[!rownames(points) %in% genomes] <- "all"
    points$compartment <- factor(points$compartment, levels=c("all", "leaf", "root", "soil"))
    points <- droplevels(points)

    colors2 <- data.frame(compartment=c("all", "leaf", "root", "soil"), 
                          cols=c("lightgrey", c_dark_green, c_red, c_red))
    c2 <- as.character(colors2$cols[match(levels(points$compartment), colors2$compartment)])
    
    p2 <- ggplot(points, aes(x=x, y=y, color=compartment, shape=taxonomy)) +
                 geom_point(alpha=.8) +
                 scale_shape_manual(values=c(16, 17, 7, 15, 3)) +
                 scale_colour_manual(values=c2, guide=F) +
                 labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
                      y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
                      title="functional distances between isolates") + 
                 main_theme +
                 theme(legend.title=element_blank(),
                        legend.text=element_text(size=8))
    
    return(arrangeGrob(p1, p2, ncol=1, heights=c(3.4, 4.9), widths=6))

}

genPangenome <- function(k) {

    # pan-genome histogram
    
    ko.all <- data.frame(genome=NULL, ko=NULL) 
    
    for (g in genomes) {
        
        ko <- read.table(paste(annotation.dir, g, ".ko", sep=""), 
                         fill=T, header=F, sep="\t",
                         col.names=c("peg", "ko"))[, 2]
        ko.genome <- data.frame(genome=g, ko=ko)
        ko.all <- rbind(ko.all, ko.genome)
    
    }
    
    ko.table <- table(ko.all)
    perc_annotated <- (sum(ko.table[, -1]) * 100) / sum(ko.table)
    
    ko.table <- ko.table[, -1]
    ko.table <- (ko.table > 0) * 1
    ko.table <- as.vector(table(colSums(ko.table)))
    
    df <- data.frame(n_genomes=1:length(ko.table), count=ko.table)
    df$count <- df$count / sum(df$count)
    df$component <- "shell"
    df$component[1] <- "singletons"
    df$component[df$n_genomes==length(genomes)] <- "core"
    
    p3 <- ggplot(df, aes(x=n_genomes, y=count, fill=component)) +
                 geom_histogram(alpha=1, size=.5, stat="identity", color="black") +
                 scale_y_continuous(labels = percent) +
                 #~ scale_x_continuous(breaks=0:length(genomes)) +
                 scale_fill_manual(values=c("black", "grey", "white")) +
                 labs(x="number of genomes", y="annotated gene families",
                      title="pan-genome histogram") +
                 main_theme +
                 theme(legend.title=element_blank())
    
    # pan-genome stacked barplot
    
    df <- aggregate(df$count, by=list(df$component), FUN=sum)
    colnames(df) <- c("component", "count")
    df$blank <-"all"
    
    p4 <- ggplot(df, aes(x=blank, y=count, fill=component)) +
                 geom_bar(alpha=1, size=.5, stat="identity", color="black", width=.4) +
                 scale_y_continuous(labels=percent) +
                 scale_fill_manual(values=c("black", "grey", "white")) +
                 labs(x="", y="annotated gene families",
                      title="pan-genome distribution") +
                 coord_flip() +
                 main_theme +
                 theme(legend.position="none",
                       axis.ticks.y=element_blank(),
                       axis.text.y=element_blank(),
                       axis.line.y=element_blank())
     
    return(arrangeGrob(p3, p4, ncol=1, heights=c(6.3, 1.7), widths=5.7))

}

genPlot <- function(k) {

    # blank place-holder for the tree and cluster name and taxonomy

    taxonomy <- tax$Family[match(genomes[1], tax$genome)]

    title <- paste("Cluster ", k, "\n(", taxonomy, ")", "\n", sep="")

    bp <- ggplot() + 
          geom_blank(aes(1, 1)) +
          labs(title=title) +
          theme(plot.background=element_blank(),
                panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.border=element_blank(),
                panel.background=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank(),
                axis.text.x=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks=element_blank())

    # arrange PDF
                 
    pg <- arrangeGrob(p.clustering, p.pangenome, ncol=2, nrow=1, heights=8.3, widths=c(6, 5.7))
    pg2 <- arrangeGrob(bp, pg, ncol=2, nrow=1, heights=8.3, widths=c(5, 11.7))
    ggsave(pg2, file=paste(figures.dir, "/amphora_tmp.pdf", sep=""), height=8.3, width=16.7) 

    system(paste("pdftk ", figures.dir, "/amphora_tmp_tree.pdf background ",
                 figures.dir, "/amphora_tmp.pdf output ",
                 figures.dir, "/cluster_", k, "_tmp.pdf", sep=""))

}

### cluster analysis plots

message("generating cluster analysis PDF...")

for (cluster in unique(groups_ids)) {

    message(paste("Cluster ", cluster, sep=""))
        
    genomes <- names(groups_ids[groups_ids==cluster])
        
    if (length(genomes) >= 5) {
        
       df <- genTree(cluster)
       p.clustering <- genDensities(cluster, df, mds)
       p.pangenome <- genPangenome(cluster) 
       genPlot(cluster) 
       
    }

}

# concatenate cluster analysis PDFs

groups_table <- table(groups_ids)
command <- gsub("$", "_tmp.pdf", names(groups_table)[groups_table >= 5])
command <- gsub("^", paste(figures.dir, "cluster_", sep=""), command)
command <- Reduce(command,f=paste)
system(paste("pdftk ", command, " output ", figures.dir, "cluster_analysis.pdf", sep=""))

# cleanup temporary files

system(paste("rm -f ", amphora.dir, "amphora_tmp*", sep=""))
system(paste("rm -f ", figures.dir, "amphora_tmp*", sep=""))
system(paste("rm -f ", figures.dir, "cluster_*_tmp.pdf", sep=""))

