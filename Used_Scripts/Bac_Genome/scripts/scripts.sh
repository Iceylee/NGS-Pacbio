#!/bin/bash
mkdir scripts && cd scripts

#contig 1
echo 'BEGIN { FS = " # " ; OFS="\t" }
/>/ {

        i = i + 1
        start = $2
        end = $3
        if ($4 == 1){
                chain = "+"
        }
        else{
                chain = "-"
        }
        printf (">%s_%06d locus=%d:%d:%s\n",name,i,start,end,chain)}
!/>/ {print $0}' > pep1.awk

echo 'BEGIN {FS="\t";OFS="\t" }
!/\#/{
        split($9,x,";")
        split(x[1],y,"_")

        printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_%08d\n",$1,$2,$3,$4,$5,$6,$7,$8,name,y[2])
}
/\#/{
        print($0)

}' > gff1.awk

#multiple contigs
echo 'BEGIN { FS = " # " ; OFS="\t" }
/>/ {

        i = i + 1
        split($1,x,"_")
        contig = "unitig_"x[2]
        start = $2
        end = $3
        if ($4 == 1){
                chain = "+"
        }
        else{
                chain = "-"
        }
        printf (">%s_%06d locus=%s:%d:%d:%s\n",name,i,contig,start,end,chain)}
!/>/ {print $0}' > pep2.awk

echo 'BEGIN {FS="\t";OFS="\t" }
!/\#/{
        split($9,x,";")
        i = i + 1

        printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=%s_%08d\n",$1,$2,$3,$4,$5,$6,$7,$8,name,i)
}
/\#/{
        print($0)

}' > gff2.awk


###
echo '#!/usr/bin/python

# coding=utf-8

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
    output_file.write(line.strip() + "\t" + anno + "\n")'>db_anno_extract.py




cd ..