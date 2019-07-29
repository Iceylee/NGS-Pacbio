#!/usr/bin/python

# coding=utf-8

#python combine2_by_key.py [inputfile1] [inputfile2][outputfile] [file1_key] [file2_key] [file2_result]
#python combine2_by_key.py genomeNum_swissprot.id swiss_GO.id gene_swiss_GO.out 2 1 2
'''
1 2 3 4 5
参数1：作为query的表
参数2：作为查询库的表
参数3：输出结果
参数4：表1的查询id
参数5：表2与参数4对应的id
参数6：表2中要查询的结果
'''

import sys,re

input_file1 = open(sys.argv[1])#"genomeNum_swissprot.id"
input_file2 = open(sys.argv[2])#"swiss_GO.id"
output_file = open(sys.argv[3],"w")#"gene_swiss_GO.out"
file1_key = int(sys.argv[4])
file2_key = int(sys.argv[5])
file2_result = int(sys.argv[6])

#test
# input_file1 = open("genomeNum_swissprot.id")#"genomeNum_swissprot.id"
# input_file2 = open("swiss_GO.id")#"swiss_GO.id"
# output_file = open("gene_swiss_GO.out","w")#"gene_swiss_GO.out"
# file1_key = 2#"swiss"
# file2_key =1#"swiss"
# file2_result =2#"GO"

a_dict = {}
for line in input_file2:
    line = line.strip()
    text = line.split("\t")
    i1 = text[file2_key-1]
    i2 = text[file2_result-1]
    a_dict[i1] = i2

    
    
    
for line in input_file1:
    line = line.strip()
    index = line.split("\t")[file1_key-1]
    if index in a_dict.keys():
        result = a_dict[index]
    else:
        result = " "
    output_file.write(line.strip() + "\t" + result + "\n")

output_file.close()