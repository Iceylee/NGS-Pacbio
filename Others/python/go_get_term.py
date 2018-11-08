'''
go列为多个GO，分号间隔；每个都需要查询得到term
gene_swiss_GO.id

'''

input_file1 = open("gene_swiss_GO.id")
input_file2 = open("go_term_class.tab")
out_file = open("gene_GOterm.out","w")

id_term_dict = {}
for line in input_file2:
    line = line.strip()
    GO_id = line.split("\t")[0]
    GO_term = line.split("\t")[1]
    id_term_dict[GO_id] = GO_term

for line in input_file1:
    line = line.strip()
    text_list = line.split("\t")
    gene_id = text_list[0]
    if len(text_list)==2:
        id_term = " "
    else:
        go_id_all = (text_list[2])
        go_id_list = go_id_all.split(";")
        id_term_list = []
        for go_id in go_id_list:
            go_id = go_id.strip()
            if go_id in id_term_dict.keys():
                go_term = id_term_dict[go_id]
            else:
                go_term = ""
            id_plus_term = go_id + "#" +go_term
            id_term_list.append(id_plus_term)
                
        id_term = (";").join(id_term_list)
    out_file.write(gene_id + "\t"+ id_term + "\n")
    
out_file.close()   