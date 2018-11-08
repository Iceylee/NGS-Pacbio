
# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list = ls())

# load libraries

library(utils, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)
library(ape, quietly=T, warn.conflicts=F)
library(Biostrings, quietly=T, warn.conflicts=F)

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

# load paths to project directories
source("paths.R")

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "wgs_taxonomy.txt", sep="")
ss.list.file <- paste(data.dir, "ss_table.txt", sep="")
photo.list.file <- paste(data.dir, "photosynthesis_table.txt", sep="")
pathways.list.file <- paste(data.dir, "bacterial_pathways.txt", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)
ss.list <- read.table(ss.list.file, sep="\t", header=T)
photo.list <- read.table(photo.list.file, sep="\t", header=T)[, 1]
bact.pathways <- read.table(pathways.list.file, sep="\t", header=T)[, 1]

### KEGG pathway abundance table

message("generating functional category tables...")

unnanotated.table <- data.frame(genome=NULL, perc=NULL)
ko.table <- data.frame(genome=NULL, pathway=NULL, perc=NULL)

threshold <- 0 / 100
        
pb <- txtProgressBar(min=1, max=length(mapping$ID), style=3)
i <- 1

for (genome in mapping$ID) {

    setTxtProgressBar(pb, i)
    i <- i + 1

    # read genome annotation
    ko.genome <- read.table(paste(annotation.dir, genome, "_KEGG.txt", sep=""), sep="\t")[, 3]
    ko.genome <- table(ko.genome)
    
    # calculate percentage of unnanotated proteins
    idx <- which(names(ko.genome)=="")
    unannotated.genome <- data.frame(genome=genome, perc=ko.genome[idx] / sum(ko.genome))
    unnanotated.table <- rbind(unnanotated.table, unannotated.genome)

    # remove unknown proteins and normalize
    ko.genome <- ko.genome[-idx]
    ko.genome <- ko.genome / sum(ko.genome)

    # keep only relevant pathways
    ko.genome <- ko.genome[names(ko.genome) %in% bact.pathways]

    # remove all pathways bellow threshold perc. abundance
    ko.genome <- ko.genome[ko.genome > threshold]

    # add genome abundances to table
    ko.genome <- data.frame(genome=genome, pathway=names(ko.genome), perc=as.vector(ko.genome))
    ko.table <- rbind(ko.table, ko.genome)

}

close(pb)

# generate overall mean functional abundance table

pathway.table <- aggregate(ko.table$perc, by=list(ko.table$pathway), FUN=sum)
colnames(pathway.table) <- c("pathway", "perc")
pathway.table$perc <- pathway.table$perc / length(unique(ko.table$genome))

### secretion systems tree figure

message("generating table of genomic features...")

kos <- data.frame(genome=NULL, KO=NULL)

pb <- txtProgressBar(min=1, max=length(mapping$ID), style=3)
i <- 1

for (genome in mapping$ID) {

    setTxtProgressBar(pb, i)
    i <- i + 1

    f <- paste(annotation.dir, genome, "_KEGG.txt", sep="")
    sec.genome <- read.table(f, sep="\t")[, 5]
    sec.genome <- table(sec.genome)
    sec.genome <- data.frame(genome=genome, ko=sec.genome)

    kos <- rbind(kos, sec.genome)

}

close(pb)

sec.table <- kos[kos[, 2] %in% ss.list$ko, 1:2]
colnames(sec.table) <- c("genome", "ko")
sec.table <- merge(sec.table, ss.list, by="ko")
sec.table <- table(sec.table[, 2:3])
sec.table <- data.frame(sec.table)

message("plotting Fig. SX (secretion systems tree figure)...")

# read tree species tree
# NOTE: manually root the tree for plotting

tree  <- read.tree(paste(amphora.dir, "amphora_rooted_names.tree", sep=""))

sec.table$presence <- sec.table$Freq > 2

# NOTE: for some bacteria (e.g. most gramm+) the TatA and TatC 
# components alone are sufficient for a functional Tat sec. system
sec.table$presence[sec.table$system=="TAT" & sec.table$Freq >= 2] <- T

# NOTE: only one gene might be sufficient for a functioning T5SS
sec.table$presence[sec.table$system=="T5SS" & sec.table$Freq >= 1] <- T

idx <- match(tree$tip.label, mapping$Strain)
sec.table$genome <- factor(sec.table$genome, levels=mapping$ID[idx])

colors <- data.frame(tax=c("Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria",
                           "Deinococcus-Thermus", "Actinobacteria", "Bacteroidetes", "Firmicutes"),
                     cols=c(c_dark_green, c_sea_green, c_green,
                            "black", c_yellow, c_blue, c_orange))

shapes <- data.frame(compartment=c("root", "soil", "leaf"), shape=c(21, 22, 24))

tax <- as.character(taxonomy$Phylum[match(tree$tip.label, taxonomy$isolate_ID)])
idx <- tax=="Proteobacteria"
tax[idx] <- as.character(taxonomy$Class[match(tree$tip.label[idx], taxonomy$isolate_ID)]) 
tips.col <- as.character(colors$cols[match(tax, colors$tax)])

comp <- mapping$compartment[match(tree$tip.label, mapping$Strain)]
tips.shape <- shapes$shape[match(comp, shapes$compartment)]

tips.cex <- rep(.6, length(tree$tip.label))

plot_btree(tree, output.file=paste(figures.dir, "/spp_tree.pdf", sep=""),
           cex.labels=.2, label.offset=.1,
   	       type="phylogram", use.edge.length=T, edge.width=1,
           align.labels=T, color.labels="black", color.alg.lines="transparent",
           tips.col=tips.col, tips.shape=tips.shape, tips.cex=tips.cex,
           margins=c(.1, 0, .1, 2), height=10, width=8)

p1 <- ggplot(sec.table, aes(x=system, y=genome, fill=presence)) +
      geom_tile() +
      scale_fill_manual(values=c("transparent", c_red)) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 20, 0, 0), "mm"))

ggsave(p1, file=paste(figures.dir, "ss.heatmap.pdf", sep=""), height=10, width=5)

### tree figure
# Carbohydrate metabolism tree figure
# Xenobiotics biodegradation and metabolism tree figure
# Photosynthesis

message("plotting Fig. SX and Fig. SX (pathway abundances tree figures)...")

pathway <- "Carbohydrate metabolism"
idx <- which(ko.table$pathway==pathway)
carb.abundances <- ko.table$perc[idx]

pathway <- "Xenobiotics biodegradation and metabolism"
idx <- which(ko.table$pathway==pathway)
xeno.abundances <- ko.table$perc[idx]

barplots.height <- data.frame(carb.abundances, xeno.abundances)

rownames(barplots.height) <- mapping$Strain[match(levels(ko.table$genome), mapping$ID)]
barplots.height <- barplots.height[match(tree$tip.label, rownames(barplots.height)), ]

annotated.proteins <- rowSums(table(kos[, 1:2]))
photo.table <- kos[kos[, 2] %in% photo.list, 1:2]
colnames(photo.table) <- c("genome", "ko")
photo.abundances <- rowSums(table(photo.table)) / annotated.proteins
names(photo.abundances) <- mapping$Strain[match(names(photo.abundances), mapping$ID)]
photo.abundances <- photo.abundances[match(tree$tip.label, names(photo.abundances))]
barplots.height$photo.abundances <- photo.abundances * 10

barplots.col <- matrix("black", nrow=dim(barplots.height)[1], ncol=dim(barplots.height)[2])
idx <- which(tree$tip.label %in% mapping$Strain[mapping$compartment=="leaf"])
barplots.col[idx, ] <- c_dark_green
idx <- which(tree$tip.label %in% mapping$Strain[mapping$compartment=="root"])
barplots.col[idx, ] <- c_red
idx <- which(tree$tip.label %in% mapping$Strain[mapping$compartment=="soil"])
barplots.col[idx, ] <- c_dark_brown

tree$tip.label <- as.character(taxonomy$Family[match(tree$tip.label, taxonomy$isolate_ID)])

plot_btree(tree, type="phylogram",
           cex.labels=.2, label.offset=.1,
           align.labels=T, color.labels="black", color.alg.lines="transparent",
           tips.col=tips.col, tips.shape=tips.shape, tips.cex=tips.cex,
           barplots=T, barplots.height=barplots.height, barplots.width=.5,
           barplots.col=barplots.col, barplots.offset=.2,
           margins=c(.1, 0, .1, 2), height=10, width=8,
           output.file=paste(figures.dir, "/tree_figure.pdf", sep=""))

