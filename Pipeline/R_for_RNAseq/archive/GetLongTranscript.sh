#!/bin/bash

if [ $# != 2 ];then
        echo
        echo "          USAGE: $0 InFasta OutFasta"
        echo "          e.g.: $0 Trinity_CD-HIT_0.9.fa Trinity_CD-HIT_0.9.long.fa"
        echo
        exit 1;
fi

InFasta=$1
OutFasta=$2

#转录本ID 长度
seqkit fx2tab -l -g -n -i -H $InFasta |cut -f 1,4 |sort -k1,1 >temp1

#新增第三列为转录本ID（去掉最后的_i1)
#按转录本名称和长度从大到小排序
#第二列相同则合并长度列。此时第一列保留为唯一最长的转录本ID
cat temp1 |awk 'BEGIN{FS=OFS="\t"}
{
	split($1,x,"_")
	print $1,x[1]"_"x[2]"_"x[3]"_"x[4],$2
}' | sort -k2,2 -k3,3nr |awk 'BEGIN { FS=OFS="\t" } 
{
    curr = $2
    if (curr == prev) {
        rec = rec ";" $3
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }' |cut -f 1|sed '1d' > id.txt

cat $InFasta|seqkit grep -f id.txt > $OutFasta

rm -f id.txt temp1