## ---- echo=FALSE---------------------------------------------------------
library(knitr)
opts_chunk$set(tidy=TRUE,message=FALSE)

## ------------------------------------------------------------------------
library(tximportData)
dir <- system.file("extdata", package="tximportData")
list.files(dir)

## ------------------------------------------------------------------------
samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
samples
files <- file.path(dir, "salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
all(file.exists(files))

## ------------------------------------------------------------------------
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype="GENEID")
df <- select(txdb, keys=k, keytype="GENEID", columns="TXNAME")
tx2gene <- df[,2:1] # tx ID, then gene ID

## ------------------------------------------------------------------------
tx2gene <- read.csv(file.path(dir, "tx2gene.csv"))
head(tx2gene)

## ------------------------------------------------------------------------
library(tximport)
library(readr)
txi <- tximport(files, type="salmon", tx2gene=tx2gene)
names(txi)
head(txi$counts)

## ------------------------------------------------------------------------
txi.tx <- tximport(files, type="salmon", txOut=TRUE, tx2gene=tx2gene)

## ------------------------------------------------------------------------
txi.sum <- summarizeToGene(txi.tx, tx2gene)
all.equal(txi$counts, txi.sum$counts)

## ------------------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
txi.salmon <- tximport(files, type="salmon", tx2gene=tx2gene)
head(txi.salmon$counts)

## ------------------------------------------------------------------------
files <- file.path(dir,"sailfish", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
txi.sailfish <- tximport(files, type="sailfish", tx2gene=tx2gene)
head(txi.sailfish$counts)

## ----eval=FALSE----------------------------------------------------------
#  txi <- tximport("quant.sf", type="none", txOut=TRUE,
#                  txIdCol="Name", abundanceCol="TPM",
#                  countsCol="NumReads", lengthCol="Length",
#                  importer=function(x) read_tsv(x, skip=8))

## ------------------------------------------------------------------------
files <- file.path(dir,"salmon_gibbs", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
txi.inf.rep <- tximport(files, type="salmon", txOut=TRUE)
names(txi.inf.rep)
names(txi.inf.rep$infReps)
dim(txi.inf.rep$infReps$sample1)

## ------------------------------------------------------------------------
files <- file.path(dir, "kallisto_boot", samples$run, "abundance.h5")
names(files) <- paste0("sample",1:6)
txi.kallisto <- tximport(files, type="kallisto", txOut=TRUE)
head(txi.kallisto$counts)

## ------------------------------------------------------------------------
names(txi.kallisto)
names(txi.kallisto$infReps)
dim(txi.kallisto$infReps$sample1)

## ------------------------------------------------------------------------
files <- file.path(dir, "kallisto", samples$run, "abundance.tsv")
names(files) <- paste0("sample",1:6)
txi.kallisto.tsv <- tximport(files, type="kallisto", tx2gene=tx2gene)
head(txi.kallisto.tsv$counts)

## ------------------------------------------------------------------------
files <- file.path(dir,"rsem", samples$run, paste0(samples$run, ".genes.results"))
names(files) <- paste0("sample",1:6)
txi.rsem <- tximport(files, type="rsem")
head(txi.rsem$counts)

## ---- results="hide", messages=FALSE-------------------------------------
library(edgeR)

## ------------------------------------------------------------------------
cts <- txi$counts
normMat <- txi$length
normMat <- normMat / exp(rowMeans(log(normMat)))
library(edgeR)
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
y <- DGEList(cts)
y$offset <- t(t(log(normMat)) + o)
# y is now ready for estimate dispersion functions
# see edgeR User's Guide

## ---- results="hide", messages=FALSE-------------------------------------
library(DESeq2)

## ------------------------------------------------------------------------
sampleTable <- data.frame(condition=factor(rep(c("A","B"),each=3)))
rownames(sampleTable) <- colnames(txi$counts)

## ------------------------------------------------------------------------
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~ condition)
# dds is now ready for DESeq()
# see DESeq2 vignette

## ------------------------------------------------------------------------
files <- file.path(dir,"salmon", samples$run, "quant.sf")
names(files) <- paste0("sample",1:6)
txi <- tximport(files, type="salmon",
                tx2gene=tx2gene,
                countsFromAbundance="lengthScaledTPM")
library(limma)
y <- DGEList(txi$counts)
y <- calcNormFactors(y)
design <- model.matrix(~ condition, data=sampleTable)
v <- voom(y, design)
# v is now ready for lmFit()
# see limma User's Guide

## ------------------------------------------------------------------------
sessionInfo()

