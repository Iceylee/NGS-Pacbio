### For human

PathView <- function(gene_exp, path_info){
  
    Gene <- read.table(gene_exp, sep="\t", header = T, row.names = 1)
  
    gene_data_matrix <- data.matrix(Gene)
  
  
    Path <- read.table(path_info, header = T, row.names = 1, sep="\t")
    
    KEGG_path_id <- (rownames(Path))
    
    cc <- function(KEGG_path_id){
        pv.out <- pathview(gene.data = gene_data_matrix[, 2], pathway.id = KEGG_path_id, 
                          species = "hsa", kegg.native = TRUE,
                          gene.idtype = "ENSEMBL")
    }
 
    sapply(KEGG_path_id, cc)  
}

sigfiles = list.files(".",pattern="sig_genes_exprData.txt")

sapply(sigfiles, PathView, path_info = "KEGG_out.txt")



## gene.idtype.list

### denovo
PathView <- function(gene_exp, path_info){
  
  Gene <- read.table(gene_exp, sep="\t", header = F)
  
  index <- duplicated(Gene$V1)
  
  gene_data_dup <- Gene[!index,]
  
  row.names(gene_data_dup) <- gene_data_dup$V1
  
  gene_data_matrix <- data.matrix(gene_data_dup)
  
  Path <- read.table(path_info, header = T, row.names = 1, sep="\t")
  
  KEGG_path_id <- (rownames(Path))
  
  View(gene_data_matrix)
  View(KEGG_path_id)
  
  cc <- function(KEGG_path_id){
    pv.out <- pathview(gene.data = gene_data_matrix[, 4], pathway.id = KEGG_path_id, 
                       species = "ko", kegg.native = TRUE)
  }
  
  sapply(KEGG_path_id, cc) 

}


sigfiles = list.files(".",pattern="LvsA_sig_pathview.txt")

sapply(sigfiles, PathView, path_info = "4.KEGG_Enrichment/LvsA.KEGG.Enrichment_out.txt")


