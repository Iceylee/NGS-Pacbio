#!/usr/bin/env python
# coding: utf-8

### Date: 2018.7.16
### Author: zss
### For: Get otu's seq 

import sys
from Bio import SeqIO


if len(sys.argv) != 3:
    print "\n\tUsage: python %s otu_table.txt rep_set.fna\n\n" % (sys.argv[0])
    sys.exit(1)

otu_table = sys.argv[1]
rep_set = sys.argv[2]

QuerySet = set()
with open(otu_table, 'r') as F:
    for line in F:
        if not line.startswith('#'):
            ID = line.strip().split('\t')[0]
            QuerySet.add(ID)

for seq_record in SeqIO.parse(rep_set, "fasta"):
    for i in QuerySet:
        if i == seq_record.id:
            print ">%s\n%s" % (seq_record.description, (seq_record.seq))
