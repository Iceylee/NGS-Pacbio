
# scripts to reproduce the analysis and figures from Bai et al., 2015
#
# originally by Ruben Garrido-Oter
# garridoo@mpipz.mpg.de

# cleanup

rm(list = ls())

# load libraries

library(utils, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(MASS, quietly=T, warn.conflicts=F)
library(gridExtra, quietly=T, warn.conflicts=F)
library(scales, quietly=T, warn.conflicts=F)

options(warn=-1)

# plotting stuff

main_theme <- theme(panel.background=element_blank(),
              panel.grid=element_blank(),
              axis.line=element_line(color="black"),
              axis.ticks=element_line(color="black"),
              axis.text=element_text(colour="black", size=10),
              legend.position="top",
              legend.background=element_blank(),
              legend.key=element_blank(),
              text=element_text(family="sans"))

alpha <- .8
c_red <-             rgb(200 / 255,   0 / 255,   0 / 255, alpha)
c_dark_green <-      rgb( 50 / 255, 200 / 255, 100 / 255, alpha)
c_dark_brown <-      rgb(101 / 255,  67 / 255,  33 / 255, alpha)

# load paths to project directories
source("paths.R")

mapping.file <- paste(data.dir, "mapping.txt", sep="")
taxonomy.file <- paste(data.dir, "wgs_taxonomy.txt", sep="")

# load data

mapping <- read.table(mapping.file, sep="\t", header=T, colClasses="character")
taxonomy <- read.table(taxonomy.file, sep="\t", header=T)

### functional profiles

message("generating matrix of functional profiles...")

ko.all <- data.frame(genome=NULL, ko=NULL) 
sizes.all <- data.frame(genome=NULL, size=NULL) 

pb <- txtProgressBar(min=1, max=length(mapping$ID), style=3)
i <- 1

for (g in mapping$ID) {
 
    setTxtProgressBar(pb, i)
    i <- i + 1
   
    ko <- read.table(paste(annotation.dir, g, ".ko", sep=""),
                     fill=T, header=F, sep="\t",
                     col.names=c("peg", "ko"))[, 2]
    ko.genome <- data.frame(genome=g, ko=ko)
    ko.all <- rbind(ko.all, ko.genome)

    size.genome <- data.frame(genome=g, size=dim(ko.genome)[1])
    sizes.all <- rbind(sizes.all, size.genome)

}

close(pb)

ko.table <- table(ko.all)
ko.table <- t(ko.table[, -1])

sizes.all$perc_annotated <- colSums(ko.table) / sizes.all$size

func <- (ko.table > 0) * 1
#~ func <- ko.table

write.table(func, file=paste(data.dir, "/functional_profiles.txt", sep=""),
            sep="\t", quote=F, col.names=T, row.names=T)

# calculate pairwise functional distances

message("calculating pairwise functional distances...")

d <- 1 - cor(func)
diag(d) <- 0

### PCoA of functional distances

message("calculating functional PCoA...")

k <- 2

pcoa <- cmdscale(d, k=k, eig=T)
points <- pcoa$points
eig <- pcoa$eig

points <- as.data.frame(points)
colnames(points) <- c("x", "y")

points$compartment <- mapping$compartment[match(rownames(points), mapping$ID)] 

taxonomy$genome <- mapping$ID[match(taxonomy$isolate_ID, mapping$Strain)]
points$taxonomy <- taxonomy$Phylum[match(rownames(points), taxonomy$genome)] 

p1 <- ggplot(points, aes(x=x, y=y, color=compartment, shape=taxonomy)) +
      geom_point(alpha=.8) +
      scale_shape_manual(values=c(16, 17, 7, 15, 3)) +
      scale_colour_manual(values=c(c_dark_green, c_red, c_dark_brown)) +
      labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
           y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep="")) + 
      main_theme +
      theme(legend.title=element_blank())

### functional distances per family

message("calculating functional distances per family...")

families <- levels(taxonomy$Family)

df <- data.frame(family=NULL, distances=NULL)

pb <- txtProgressBar(min=1, max=length(families), style=3)
i <- 1

for (f in families) {

    setTxtProgressBar(pb, i)
    i <- i + 1

    idx <- which(taxonomy$Family[match(mapping$Strain[match(rownames(d), mapping$ID)],
                                       taxonomy$isolate_ID)] == f)

    distances <- d[idx]
    phylum <- taxonomy$Phylum[idx]
    compartment <- mapping$compartment[match(taxonomy$isolate_ID[idx], mapping$Strain)]
    perc <- sizes.all$perc_annotated[idx]

    if (length(idx) >= 5) {
        
        df <- rbind(df, data.frame(family=f, distances=distances,
                                   phylum=phylum, compartment=compartment,
                                   perc_annotated=perc))

    }

}

close(pb)

message("plotting Fig. 3...")

medians <- aggregate(df$distances, by=list(df$family), FUN=median)
order <- medians[sort(medians[, 2], index.return=T, decreasing=F)$ix, 1]
df$family <- factor(df$family, levels=order)

p2 <- ggplot(df, aes(x=family, y=distances, shape=phylum, color=compartment)) +
                 geom_jitter(position=position_jitter(0.35), size=.7, alpha=1) +
                 geom_boxplot(alpha=1, outlier.size=.5, size=.3, color="black", fill=NA) +
                 scale_shape_manual(values=c(16, 17, 7, 3, 3), guide=F) +
                 scale_colour_manual(values=c(c_dark_green, c_red, c_dark_brown), guide=F) +
                 labs(x="", y="functional distance") +
                 coord_flip() +
                 scale_y_continuous(limits=c(.25, .75)) +
                 main_theme +
                 theme(legend.title=element_blank(),
                       title=element_blank(),
                       axis.text.y = element_text(size=6))

p3 <- ggplot(df, aes(x=distances)) +
             geom_histogram(size=.5, alpha=1, color="black", fill="grey", binwidth=.01) +
             scale_x_continuous(limits=c(.24, .75)) +
             labs(title="pairwise functional distances between isolates of the same family", x="", y="") +
             main_theme +
             theme(legend.position="none", axis.text.y=element_text(size=8),
                   title=element_text(size=6))

p.fam <- arrangeGrob(p3, p2, ncol=1, nrow=2, heights=c(2, 6.27), widths=c(2, 5))

pg <- arrangeGrob(p1, p.fam, ncol=2, nrow=1, heights=8.27, widths=c(8, 4))
ggsave(file=paste(figures.dir, "/functional_MDS.pdf", sep=""), pg, height=8.27, width=13)

p4 <- ggplot(df, aes(x=family, y=perc_annotated, shape=phylum, color=compartment)) +
                 geom_jitter(position=position_jitter(0.35), size=.7, alpha=1) +
                 geom_boxplot(alpha=1, outlier.size=.5, size=.3, color="black", fill=NA) +
                 scale_shape_manual(values=c(16, 17, 7, 3, 3), guide=F) +
                 scale_colour_manual(values=c(c_dark_green, c_red, c_dark_brown), guide=F) +
                 labs(x="", y="percentage of annotated proteins") +
                 scale_y_continuous(labels=percent) +
                 coord_flip() +
                 main_theme +
                 theme(legend.title=element_blank(),
                       title=element_blank(),
                       axis.text.y = element_text(size=6))

p.fam <- arrangeGrob(p2, p4, ncol=2, nrow=1, heights=6.27, widths=4)
ggsave(file=paste(figures.dir, "/func_dist_v_perc_annotated.pdf", sep=""), p.fam, height=6.27, width=8)

