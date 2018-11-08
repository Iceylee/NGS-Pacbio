#!/bin/bash

if [ $# != 2 ];then
        echo 
        echo "          USAGE: $0 inputfile tabfile"
        echo "          (Note: In tabfile 1st col from ID;2nd col to ID"
        echo "          e.g.: $0 1GO_BP_out.txt ID.list"
        echo 
        exit 1;
fi


awk 'BEGIN{FS=OFS="\t"}{split($8,x,"/");for (i in x)print $1,$2,$3,$4,$5,$6,$7,x[i],$9}' $1>temp1


#对应ID
awk '
NR==FNR{FS=",";a[$1]=$2;next}
{FS=OFS="\t"}
{print $1,$2,$3,$4,$5,$6,$7,a[$8],$9}' $2 temp1  > temp2

#合并（第1列相同，合并第八列）
echo 'BEGIN { FS=OFS="\t" }
FNR==1{print $0}
FNR>1{
    curr = $1
    if (curr == prev) {
        rec = rec "/"$8
    }
    else {
        if (rec) print rec,last_nine
        rec = $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8
    }
    prev = curr
    last_nine = $9
}
END { if (rec) print rec }' > merge.awk

awk -f merge.awk temp2 > temp3

#mv temp3 $1