
'''
将clusterprofiler的三个输出文件，转换成id-GO并合并。

当前文件夹有三个输入文件 GO_BP_out.txt GO_MF_out.txt GO_CC_out.txt
输出文件为 id_GO_all.out
---
entrezID GO编号+term GO分类MF/BP/CC
---

python GOid_to_idGO.py

'''

GO_class = ["BP","CC","MF"]

out_file_name = "id_GO_all.out"
output_file=open(out_file_name,"w")

for i in GO_class:
    in_file_name = "GO_"+i+"_out.txt"
    input_file=open(in_file_name)   

    id_go_dic = {}
    for line in input_file:
        GO_num = line.split("\t")[0]
        GO_description = line.split("\t")[1]
        GO_term = GO_num + " " + GO_description
        ids = (line.strip()).split("\t")[7]
        id_list = ids.split("/")
        for id_name in id_list:
            if id_name not in id_go_dic:
                id_go_dic[id_name]=[]
            id_go_dic[id_name].append(GO_term)


    for key in id_go_dic:
        go_list = (";").join(id_go_dic[key])
        output_file.write(key + "\t" + go_list + "\t" + i + "\n")

output_file.close()