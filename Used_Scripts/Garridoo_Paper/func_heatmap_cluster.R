
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

# input data files

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "wgs_taxonomy.txt", sep="")
pathways.list.file <- paste(data.dir, "bacterial_pathways.txt", sep="")
clusters.file <- paste(amphora.dir, "clusters.txt", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)
clusters <- read.table(clusters.file, sep="\t", header=T)
bact.pathways <- read.table(pathways.list.file, sep="\t", header=F)[, 1]

### pick representatives of each cluster and generate the tree for the figure

message("generating tree of cluster representative sequences...")

reps <- data.frame(cluster=NULL, representative=NULL, root.count=NULL,
                   soil.count=NULL, leaf.count=NULL, family=NULL, phylum=NULL)

for (cluster in unique(clusters$cluster)) {

    # take first member of each cluster as a representative
    representative <- clusters$ID[clusters$cluster==cluster][1]
    rep.c <- data.frame(cluster, representative)

    genomes <- clusters$ID[clusters$cluster==cluster]

    counts <- table(mapping$compartment[match(genomes, mapping$Strain)])
    rep.c$leaf.count <- counts["leaf"]
    rep.c$root.count <- counts["root"]
    rep.c$soil.count <- counts["soil"]
    rep.c[is.na(rep.c)] <- 0

    idx <- match(genomes, taxonomy$isolate_ID)
    rep.c$family <- Reduce(x=unique(taxonomy$Family[idx]), paste)
    rep.c$phylum <- Reduce(x=unique(taxonomy$Phylum[idx]), paste)

    reps <- rbind(reps, rep.c)

}

# save table of cluster representatives
write.table(reps, file=paste(amphora.dir, "cluster_reps.txt", sep=""),
            sep="\t", col.names=T, row.names=F, quote=F)

# read amphora MSA
amphora.msa <- readDNAStringSet(paste(amphora.dir, "amphora.msa", sep=""))

# subset sequences from the amphora MSA

rep.genomes <- mapping$ID[match(reps$representative, mapping$Strain)]
seqs <- amphora.msa[names(amphora.msa) %in% rep.genomes]
names(seqs) <- mapping$Strain[match(names(seqs), mapping$ID)]

writeXStringSet(seqs, paste(amphora.dir, "cluster_reps.msa", sep=""))

# build sub-tree

system(paste("rm -f ", amphora.dir, "cluster_reps.tree", sep=""))

system(paste("FastTree -nt -gtrz ", amphora.dir, "cluster_reps.msa >> ",
             amphora.dir, "cluster_reps.tree 2>> /dev/null", sep=""))

# plot tree of cluster representative sequences
# NOTE: root the tree manually to display in the figure

tree  <- read.tree(paste(amphora.dir, "cluster_reps_rooted.tree", sep=""))

reps <- reps[match(tree$tip.label, reps$representative), ]

idx <- reps$phylum=="Proteobacteria"
reps$phylum <- as.character(reps$phylum)
reps$phylum[idx] <- as.character(taxonomy$Class[match(reps$representative[idx], taxonomy$isolate_ID)])

colors <- data.frame(tax=c("Gammaproteobacteria", "Betaproteobacteria", "Alphaproteobacteria",
                           "Deinococcus-Thermus", "Actinobacteria", "Bacteroidetes", "Firmicutes"),
                     cols=c(c_dark_green, c_sea_green, c_green,
                            "black", c_yellow, c_blue, c_orange))

tips.col <- as.character(colors$cols[match(reps$phylum, colors$tax)])

tips.shape <- rep(21, length(tree$tip.label))

tips.cex <- rep(1, length(tree$tip.label))

tree$tip.label <- apply(reps[, c(1, 3, 4, 5, 6)], 1, function(i) Reduce(x=i, paste))

plot_btree(tree, output.file=paste(figures.dir, "/cluster_tree.pdf", sep=""),
           cex.labels=.5, label.offset=.1,
   	       type="phylogram", use.edge.length=T, edge.width=2,
           align.labels=T, color.labels="black", color.alg.lines="black",
           tips.col=tips.col, tips.shape=tips.shape, tips.cex=tips.cex,
           margins=c(.1, 0, .1, 7.5), height=6.5, width=11.7)

### KEGG pathway abundance table and heatmap
# NOTE: functional categories here correspond to the 
# sub-class KEGG category, instead of proper KEGG pathways
# (which are one level down in the BRITE hierarchy)

message("generating functional abundance tables...")

unnanotated.table <- data.frame(genome=NULL, perc=NULL)
ko.table <- data.frame(genome=NULL, pathway=NULL, perc=NULL)

threshold <- 0 / 100

pb <- txtProgressBar(min=1, max=length(mapping$ID), style=3)
i <- 1

for (genome in mapping$ID) {

    setTxtProgressBar(pb, i)
    i <- i + 1

    # read genome annotation
    ko.genome <- read.table(paste(annotation.dir, genome, "_KEGG.txt", sep=""),
                            fill=T, sep="\t", header=F,
                            col.names=c("peg", "class", "subclass", "pathway", "ko", "description"))[, 3]
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

# aggregate KO abundances for each cluster

func.table <- data.frame(pathway=NULL, perc=NULL, cluster=NULL)

for (cluster in unique(clusters$cluster)) {

    genomes <- clusters$ID[clusters$cluster==cluster]
    genomes <- mapping$ID[match(genomes, mapping$Strain)]
    cluster.table <- ko.table[ko.table$genome %in% genomes, -1]
    cluster.table <- aggregate(cluster.table$perc, by=list(cluster.table$pathway), FUN=sum)
    colnames(cluster.table) <- c("pathway", "perc")
    cluster.table$perc <- cluster.table$perc / length(genomes)
    cluster.table$cluster <- cluster

    func.table <- rbind(func.table, cluster.table)

}

# order cluster levels to match the tree tips

func.table$cluster <- factor(func.table$cluster, levels=reps$cluster)

# order functions by decreasing mean rel. abundance

idx <- sort(pathway.table$perc, decreasing=T, index.return=T)$ix
ko.table$pathway <- factor(ko.table$pathway, levels=pathway.table$pathway[idx])
pathway.table$pathway <- factor(pathway.table$pathway, levels=pathway.table$pathway[idx])
func.table$pathway <- factor(func.table$pathway, levels=pathway.table$pathway[idx])

### plot heatmap

colors <- c("#030303", "#454409", "#7b7807", "#b0a609",
            "#ddd300", "#ffff04", "#ffcd00", "#ff9c0a",
            "#fd6e10", "#fe2f1d")
b <- c(-Inf, 0.00, 0.01, 0.02, 0.04, 0.06, 0.08, 0.09, 0.1, 0.12, 0.15, 1) * 100
func.table$perc_bins <- cut(func.table$perc * 100, breaks=b, right=F)

p1 <- ggplot(func.table, aes(x=pathway, y=cluster, fill=perc_bins)) +
      geom_tile() +
      scale_fill_manual(values=colors) +
      labs(x="", y="", title="") + 
      main_theme +
      theme(legend.title=element_blank(),
            legend.position="left",
            axis.text.x=element_text(angle=45, hjust=1, size=10),
            axis.text.y=element_blank(),
            axis.ticks=element_blank(),
            axis.line=element_blank(),
            plot.margin=unit(c(0, 0, 0, 0), "mm"))

### function enrichment tests (root v. leaf .v soil) and boxplot figure

ko.table$compartment <- mapping$compartment[match(ko.table$genome, mapping$ID)]
ko.table$associated <- mapping$associated[match(ko.table$genome, mapping$ID)]

ko.table.root <- ko.table[ko.table$compartment=="root", ]
mean <- aggregate(ko.table.root$perc, by=list(ko.table.root$pathway), FUN=mean)
sd <- aggregate(ko.table.root$perc, by=list(ko.table.root$pathway), FUN=sd)
ko.table.root <- data.frame(pathway=mean[, 1], mean=mean[, 2], sd=sd[, 2], compartment="root")

ko.table.soil <- ko.table[ko.table$compartment=="soil", ]
mean <- aggregate(ko.table.soil$perc, by=list(ko.table.soil$pathway), FUN=mean)
sd <- aggregate(ko.table.soil$perc, by=list(ko.table.soil$pathway), FUN=sd)
ko.table.soil <- data.frame(pathway=mean[, 1], mean=mean[, 2], sd=sd[, 2], compartment="soil")

ko.table.leaf <- ko.table[ko.table$compartment=="leaf", ]
mean <- aggregate(ko.table.leaf$perc, by=list(ko.table.leaf$pathway), FUN=mean)
sd <- aggregate(ko.table.leaf$perc, by=list(ko.table.leaf$pathway), FUN=sd)
ko.table.leaf <- data.frame(pathway=mean[, 1], mean=mean[, 2], sd=sd[, 2], compartment="leaf")

ko.compartment.table <- rbind(ko.table.root, ko.table.soil, ko.table.leaf)

p3 <- ggplot(ko.compartment.table, aes(x=pathway, y=mean, color=compartment)) +
                 geom_pointrange(aes(ymin=mean-sd, ymax=mean+sd), size=.9,
                                 position=position_dodge(width=0.70)) +
                 scale_color_manual(values=c(c_red, c_dark_brown, c_dark_green)) +
                 scale_y_continuous(labels=percent) +
                 labs(x="", y="") +
                 main_theme +
                 theme(legend.title=element_blank(),
                       legend.position="left",
                       axis.text.x=element_blank(),
                       axis.text.y=element_text(size=9),
                       title=element_blank(),
                       plot.margin=unit(c(5, 0, 0, 3.5), "mm"))

### Mann-Whitney tests of differencial functional abundances

message("performing enrichment tests...")

p.vals <- data.frame(pathway=NULL, p.rvs=NULL, p.rvl=NULL, ci.rvs=NULL, ci.rvl=NULL)

for (pathway in bact.pathways) {

    # get vectors of relative abundances by compartment

    pathway.root <- ko.table$perc[ko.table$pathway==pathway & ko.table$compartment=="root"]
    pathway.soil <- ko.table$perc[ko.table$pathway==pathway & ko.table$compartment=="soil"]
    pathway.leaf <- ko.table$perc[ko.table$pathway==pathway & ko.table$compartment=="leaf"]

    # perform test between root and soil and root v. leaf
    rvs <- wilcox.test(pathway.root, pathway.soil, alternative="two.sided", correct=F, conf.int=T)
    rvl <- wilcox.test(pathway.root, pathway.leaf, alternative="two.sided", correct=F, conf.int=T)

    pathway.p <- data.frame(pathway=pathway, 
                            p.rvs=rvs$p.value,
                            p.rvl=rvl$p.value,
                            ci.rvs=rvs$estimate,
                            ci.rvl=rvl$estimate)

    p.vals <- rbind(p.vals, pathway.p)

}

row.names(p.vals) <- NULL

# Bonferroni multiple testing correction
# (number of functional categories * 2 two-sided groups of tests)

p.vals$p.rvs <- p.vals$p.rvs * length(bact.pathways) * 4
p.vals$p.rvl <- p.vals$p.rvl * length(bact.pathways) * 4

p.vals <- p.vals[match(levels(pathway.table$pathway), p.vals$pathway), ]

# set dignificant threshold to 0.005
alpha <- 0.05

# consider only as sig. shift in medians larger than 0.005 % rel. abundances
ci.threshold <- .5 / 100

# print results

message("significant differences between root and soil:")
print(p.vals[p.vals$p.rvs < alpha & abs(p.vals$ci.rvs) > ci.threshold, c("pathway", "p.rvs", "ci.rvs")])

message("significant differences between root and leaf:")
print(p.vals[p.vals$p.rvl < alpha & abs(p.vals$ci.rvl) > ci.threshold, c("pathway", "p.rvl", "ci.rvl")])

### pangenome analysis and stacked barplot figure

message("calculating pan-genome distribution per cluster...")

# for each cluster, calculate pangenome distribution
# (core, shell, singletons)
# with respect of percentage of annotated proteins

pangenome.table.all <- data.frame(component=NULL, perc=NULL, cluster=NULL)

for (cluster in unique(clusters$cluster)) {

    # get genome members

    genomes <- clusters$ID[clusters$cluster==cluster]
    genomes <- mapping$ID[match(genomes, mapping$Strain)]

    # if the cluster has at least than 2 members

    if (length(genomes) > 1) {

        # get cluster KO abundances

        ko.all <- data.frame(genome=NULL, ko=NULL) 
        
        for (g in genomes) {
            
            ko <- read.table(paste(annotation.dir, g, ".ko", sep=""),
                             fill=T, header=F, sep="\t",
                             col.names=c("peg", "ko"))[, 2]
            ko.genome <- data.frame(genome=g, ko=ko)
            ko.all <- rbind(ko.all, ko.genome)
        
        }
        
        ko.all <- table(ko.all)
        perc_annotated <- (sum(ko.all[, -1]) * 100) / sum(ko.all)
        
        # calculate pangenome distribution

        ko.all <- ko.all[, -1]
        ko.all <- (ko.all > 0) * 1
        ko.all <- as.vector(table(colSums(ko.all)))
        
        pangenome.table <- data.frame(n_genomes=1:length(ko.all), count=ko.all)
        pangenome.table$count <- pangenome.table$count / sum(pangenome.table$count)
        pangenome.table$component <- "shell"
        pangenome.table$component[1] <- "singletons"
        pangenome.table$component[pangenome.table$n_genomes==length(genomes)] <- "core"

        pangenome.table <- aggregate(pangenome.table$count,
                                     by=list(pangenome.table$component),
                                     FUN=sum)
        colnames(pangenome.table) <- c("component", "perc")
        pangenome.table$cluster <- cluster

        pangenome.table.all <- rbind(pangenome.table.all, pangenome.table)

    } else {

        # if cluster is a singleton, set core and shell to 0
        
        pangenome.table <- data.frame(component=c("singletons", "shell", "core"),
                                      perc=c(0, 0, 1),
                                      cluster=cluster)
        pangenome.table$cluster <- cluster
        pangenome.table$perc[1:2] <- 0
        pangenome.table$perc[3] <- 1
        
        pangenome.table.all <- rbind(pangenome.table.all, pangenome.table)

    }

}

# plotting

message("plotting Fig. 4...")

pangenome.table.all$cluster <- factor(pangenome.table.all$cluster, levels=reps$cluster)

p4 <- ggplot(pangenome.table.all, aes(x=cluster, y=perc, fill=component)) +
             geom_bar(alpha=1, size=.5, stat="identity", color="black", width=.7) +
             scale_y_continuous(labels=percent) +
             scale_fill_manual(values=c("black", "grey", "white")) +
             labs(x="", y="", title="") +
             coord_flip() +
             main_theme +
             theme(legend.position="top",
                   legend.title=element_blank(),
                   axis.title.y=element_blank(),
                   axis.ticks.y=element_blank(),
                   axis.line.y=element_blank(),
                   plot.margin=unit(c(80, 5, 47, 0), "mm"))
     
pg1 <- arrangeGrob(p3, p1, ncol=1, nrow=2, heights=c(3.5, 8.5), widths=12)
pg2 <- arrangeGrob(pg1, p4, ncol=2, nrow=1, heights=13, widths=c(12, 3))

ggsave(file=paste(figures.dir, "/heatmap_function_cluster.pdf", sep=""), pg2, height=13, width=15)

