#!/usr/bin/python

# coding=utf-8

#python extractPoly.py [mafFile] [outflile] [ref-name]
#python extractPoly.py ref_qry.maf out.out NZ_CP014343.1

import sys, re

maf = open(sys.argv[1],'rU')
out = open(sys.argv[2],'w')
ref_name = sys.argv[3]

for line in maf:
    if re.search(r"^s", line) and not re.search(ref_name, line):
            acc = line.rstrip().split()[1]
            break
maf.seek(0,0)

###

head = ("Position\t"+ref_name+"\t%s\n")%acc

out.write(head)
tagNum = 0
tag = 0
tag2 = 0
ruf = ""
seq = ""
seqTag = []
for line in maf:
    if re.search(r"^a",line):
            tag = 0
            tag2 = 0
    elif re.search(ref_name,line):
            tag = 1
            _line = line.rstrip().split()
            start_ruf = int(_line[2])
            end_ruf = int(_line[3])+1-start_ruf
            ruf = _line[6]

    elif re.search(r"^s",line) and not re.search(ref_name,line):
            tag2 = 1
            _line = line.rstrip().split()
            seq = _line[6]
    else:
            tag = 3
    if tag == 1 and tag2 == 1:
        seqTag = []  #清空？？
        for i in range(len(ruf)):
            if ruf[i]!=seq[i]:
                seqTag.append(ruf[i])
                count = seqTag.count("-")
                out.write(str(start_ruf+i-count)+"\t"+ruf[i]+"\t"+seq[i]+"\n")
                #tagNum += 1
                #print (tagNum)
maf.close()
out.close()