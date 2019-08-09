#!/usr/bin/R

args<-commandArgs(T)

if (length(args)!=1) {
	print ("")
	print("Usage: Rscript /data/script/Get_OrgDb_KEGG_Name.R \"Homo sapiens\"")
	print("")
  stop(call.=FALSE)
} 

SpeciesFullName <- args[1]

print (paste("Searching for species:",SpeciesFullName,"...",sep=""))

suppressMessages(library(clusterProfiler))
suppressMessages(library(DOSE))
suppressMessages(require(AnnotationHub))


hub <- AnnotationHub()

org <- hub[hub$rdataclass == "OrgDb",]
hm <- query(org, SpeciesFullName)
SqlitDB <- hm$title
#SpeciesName <- Query[Query$title == SqlitDB]$species

Ah_id <- (hub[hub$title == SqlitDB])$ah_id

KEGG_id <- search_kegg_organism(SpeciesFullName, by='scientific_name')

print("===================")
Ah_id
print("===================")
KEGG_id


