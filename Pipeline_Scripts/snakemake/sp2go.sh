#!/bin/bash

#从go注释结果统计go，用于作图
# sp2go.sh sp.out go_term_class.tab sp_go_merge.gaf go.anno.out go_count.out ./
out=$1
tab=$2
gaf=$3
outdir=$6


#从blast结果提取gene id 和swiss prot id
cat $out |awk 'BEGIN{FS=OFS="\t"}{split($2,x,"|");print $1,x[2]}' >$outdir/gene_swissprot.id

awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if($2 in a) print $0,a[$2]}' $gaf $outdir/gene_swissprot.id |awk 'BEGIN{FS=OFS="\t"}{split($3,x,";");for(i in x)print $1,$2,x[i]}' > $outdir/go.tab

#按照GO编号加term和class (go_term_class.tab来自go.obo，goobo_extract_info.py)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{print $0, a[$3]}' $tab $outdir/go.tab >$outdir/go_term.tab


#得到GO计数。用于作图
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{print $0, a[$3]}' $tab $outdir/go.tab | awk 'BEGIN{FS=OFS="\t"}{term[$3]=$5;class[$3]=$4;count[$3]+=1}END{for(id in count)print id,count[id],term[id],class[id]}'|sort -k2nr,2 > $5 


#合并go和term。（一个gene ID可能对应多列）
cat $outdir/go_term.tab | awk 'BEGIN{FS=OFS="\t"}{plus=$3"("$4")";print $1,plus}' >$outdir/go_plus_term.tab


#1列相同的，合并第二列（分号间隔）
echo 'BEGIN { FS=OFS="\t" }
{
    curr = $1
    if (curr == prev) {
        rec = rec ";" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }' > $outdir/merge.awk

awk -f $outdir/merge.awk $outdir/go_plus_term.tab > $4
#得到最终注释结果

rm $outdir/gene_swissprot.id $outdir/go.tab $outdir/go_term.tab $outdir/go_plus_term.tab $outdir/merge.awk




