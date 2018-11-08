
#!/usr/bin/python

# coding=utf-8

'''
python db_anno_extract.py [gene-anno-file] [blast-out-file][output-name]

去掉ID的>存为i1,将后面所有信息合并存为i2
将besthit的第二列index与i1匹配。得到i2注释作为最后一列anno添加。
'''


import sys,re

input_file1 = open(sys.argv[1])#"uniprot_sprot.anno"
input_file2 = open(sys.argv[2])#"swissprot.bestHits"
output_file = open(sys.argv[3],"w")#"swissprot.anno.out"

pro_anno = {}
for line in input_file1:
    i1 = line.split(" ")[0][1:]
    i_list = (line.strip()).split(" ")[1:]
    i2 = (" ").join(i_list)
    pro_anno[i1] = i2
    
    
for line in input_file2:
    index = line.split()[1]
    anno = pro_anno[index]
    output_file.write(line.strip() + "\t" + anno + "\n")
    
    
    