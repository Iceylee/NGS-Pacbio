## Scripts for updating specData, speciesMap and validTaxId

## Download and unpack mapping file:
## ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

## Generates specData
.processTaxNamesFile <- function(filesDir=getwd()){
##    species <- read.delim('names.dmp',header = FALSE,sep = "|")
    dest  <- file.path(filesDir, "names.dmp")
    data <-  read.delim(dest, header=FALSE, sep="\t", quote="",
                        stringsAsFactors=FALSE)
    species <- data[,seq(1, dim(data)[2], by=2)] ## Throw away 'pipe columns'
    colnames(species) <- c('tax_id','name_txt','unique_name','name_class')
    ## keep only some cols
    species <- species[,c(1:2,4)]
    ## throw away tabs from second col
    species[[2]] <- gsub('\t','',species[[2]])
    ## And the third col
    species[[3]] <- gsub('\t','',species[[3]])
    ## throw away rows where the third column doesn't say 'scientific name'
    keep <- grepl('scientific name', species[[3]])
    keep2 <- grepl('synonym', species[[3]])
    species <- species[(keep | keep2), 1:2]

    ## split second column by first space:
    rawSpec <- species[[2]]
    spltSpec <- strsplit(rawSpec, split=" ")
    genusDat <- unlist(lapply(spltSpec, function(x){x[1]}))
    .getRest <- function(x){
        if(length(x) > 1){
            return(paste(x[2:length(x)], collapse=" "))
        }else{
            return(NA)
        }
    }
    speciesDat <- unlist(lapply(spltSpec, .getRest))
    specData <- data.frame(tax_id=as.integer(species[[1]]), ## integer
                           genus=as.factor(genusDat),       ## factor
                           species=speciesDat,              ## character
                           stringsAsFactors=FALSE)
    save(specData, file='specData.rda', compress="xz")
}

## Generates speciesMap and validTaxIds
.processSpeciesMapData <- function(){
    con <- file('names.dmp')
    species <- readLines(con)
    close(con)
    splt <- strsplit(species, split='\\t\\|\\t')
    ## Throw away elements where column 4 is not 'scientific name' or 'synonym'
    idx1 <- unlist(lapply(splt, function(x){grepl('scientific name', x[4])}))
    idx2 <- unlist(lapply(splt, function(x){grepl('synonym', x[4])}))
    idx <- idx1 | idx2
    splt <- splt[idx]
    ## and keep only 1st two elements
    taxon <-  as.integer(unlist(lapply(splt, function(x){x[1]})))
    species <- unlist(lapply(splt, function(x){x[2]}))
    speciesMap <- data.frame(taxon,    ## integer
                             species,  ## character
                             stringsAsFactors=FALSE)
    save(speciesMap, file='speciesMap.rda', compress="xz")

    ## Then get the valid Tax IDs.
    validTaxIds <- unique(speciesMap$taxon)  ## integer
    save(validTaxIds, file='validTaxIds.rda', compress="xz")
}
