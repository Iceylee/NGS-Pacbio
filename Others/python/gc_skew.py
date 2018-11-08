#!/usr/bin/python

# coding=utf-8

'''
gc_skew.py [input-fasta][output-txt]

'''


import sys,re

input_file_name = sys.argv[1]#"polish_assembly.fasta"
output_file = open(sys.argv[2],"w")#"GC_skew.txt"

from Bio import SeqIO, SeqUtils

rec = SeqIO.read(input_file_name, 'fasta')
gc = SeqUtils.GC_skew(rec.seq, 1000)
start = 1
for gc_value in gc:
    end = start + 1000
    if end > len(rec):
        end = len(rec)
    output_file.write(str(start)+"\t"+str(end)+"\t"+str(gc_value) + "\n")
    start = start + 1000
    
output_file.close()