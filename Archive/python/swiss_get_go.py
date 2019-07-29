'''
输入文件：geneID-swissID goa_uniprot数据库（包含uniprotID，go编号）
输出文件：geneID-swissID-go

根据swissID查询数据库得到GO编号

'''
input_file1 = open("genomeNum_swissprot.id")
input_file2 = open("goa_uniprot_all.gaf")
output_file = open("swissprot_go.id","w")
 
goa = {}
for line in input_file2:
    if line[0] != "!":
        
        line = line.strip()
        v1 = line.split("\t")[1]
        v2 = line.split("\t")[4]

        if v1 not in goa.keys():
            goa[v1]=[]
        goa[v1].append(v2)

    
for line in input_file1:
    line = line.strip()
    gene_id = line.split(" ")[0]
    sp_id = line.split(" ")[1]
    if sp_id in goa.keys():
        go = goa[sp_id]#go is a list    
    else:
        go = []
    go_all = (";").join(go)
    output_file.write(gene_id + "\t" + sp_id + "\t" + go_all + "\n")
    
    
output_file.close()