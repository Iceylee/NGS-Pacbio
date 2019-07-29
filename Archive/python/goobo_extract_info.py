#下载的go.obo
#得到 go+term+分类 的表格

input_file = open("go.obo")
out_file = open("go_term_class.tab","w")

list_go = []
flag = False
for line in input_file:
    line = line.strip()
    if line[0:6] == "[Term]":
        flag = True
    if line[0:2] == "id":
        list_go.append(line[4:])
    if line[0:5] == "name:":
        list_go.append(line[6:])
    if line[0:9] == "namespace":
        list_go.append(line[11:])
        if flag:
            out_file.write(list_go[0]+ "\t" + list_go[1] + "\t" + list_go[2] + "\n")
        list_go = []
        flag = False

out_file.close()