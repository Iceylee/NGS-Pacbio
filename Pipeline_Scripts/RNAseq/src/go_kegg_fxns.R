suppressMessages(library(DOSE))
suppressMessages(library(clusterProfiler))
suppressMessages(require(AnnotationHub))
suppressMessages(library(tidyverse))

##' get KEGG number in AnnotationHub for later use in clusterProfiler.
##'
##'
##' @title getKEGGName
##' @param species a character of species full name
##' @return a vector of Ah_id and KEGG_id
##' @example getKEGGName("Homo sapiens")
##' @author Yubing Li




getKEGGName <- function (species){
    
    hub <- AnnotationHub()
    
    org <- hub[hub$rdataclass == "OrgDb",]
    hm <- query(org, species)
    SqlitDB <- hm$title
    #SpeciesName <- Query[Query$title == SqlitDB]$species
    
    Ah_id <- (hub[hub$title == SqlitDB])$ah_id
    
    KEGG_id <- search_kegg_organism(species, by='scientific_name')
    
    results <- c(Ah_id = Ah_id, KEGG_id = KEGG_id)
    return(results)
}




