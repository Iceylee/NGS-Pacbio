#作图用
#按照GO，将gene ID计数。并加上class 列
input_file = open("gene_swiss_GO.id")
output_file = open("go_gene_sum.count","w")

go_gene_dict = {}
for line in input_file:
    text_list = (line.strip()).split("\t")
    gene_id = text_list[0]
    if len(text_list) != 2:
        go_all = text_list[2]
        go_list = go_all.split(";")
        for go_id in go_list:
            if go_id not in go_gene_dict.keys():
                go_gene_dict[go_id] = []
            go_gene_dict[go_id].append(gene_id)



for go_id in go_gene_dict.keys():
    gene_count = len(go_gene_dict[go_id])
    gene_all = (";").join(go_gene_dict[go_id])
    go_id = go_id.strip()
    output_file.write(go_id + "\t" + str(gene_count) + "\t" + gene_all + "\n")
            
output_file.close()