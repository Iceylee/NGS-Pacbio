
#2.GO prep
library(AnnotationHub)
hub <- AnnotationHub()

#   1. model animal
# library(org.Hs.eg.db)
# db <- org.Hs.eg.db
#   2. download

org <- ah[ah$rdataclass == "OrgDb",]
hm <- query(org, "Homo sapiens")

query(hub, "Verticillium dahliae")
db <- hub[['AH59553']]

columns(db)
head(keys(db,keytype = "SYMBOL"))

#3.KEGG prep
library(clusterProfiler)
library(DOSE)
# kegg_code      scientific_name          common_name
search_kegg_organism('human', by='common_name')

#