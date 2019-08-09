cd GO
#ln -s ~/analysis/database/goa_uniprot_all_duplicates_moved.gaf ./gaf
ln -s ../bestHits/sp.bestHits ./
ln -s ~/analysis/database/go_term_class.tab ./
ln -s ~/analysis/database/sp_go_merge.gaf ./


#从blast结果提取gene id 和swiss prot id
cat sp.bestHits |awk 'BEGIN{FS=OFS="\t"}{split($2,x,"|");print $1,x[2]}' >gene_swissprot.id

#按照swissprot ID 对应得到go编号。输出为geneID+gaf所有信息
#gene ID可能对应多个swissprotID(怀疑错误 :多个geneID对应同一个spID) 因此映射应从gene到sp
#if($2 in a) 有go再输出；否则每个geneID都会输出，没有go的输出空
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{if($2 in a) print $0,a[$2]}' sp_go_merge.gaf gene_swissprot.id > go.temp #6546 （小于sp_gene.id的6817）

#拆分 每个GO一行
cat go.temp|awk 'BEGIN{FS=OFS="\t"}{split($3,x,";");for(i in x)print $1,$2,x[i]}' > go.tab

#按照GO编号加term和class (go_term_class.tab来自go.obo，goobo_extract_info.py)
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{print $0, a[$3]}' go_term_class.tab go.tab >go_term.tab

#查看重复
#cat go_term.tab|cut -f 2|sort|uniq|wc -l

#得到GO计数。用于作图
cat go_term.tab | awk 'BEGIN{FS=OFS="\t"}{term[$3]=$5;class[$3]=$4;count[$3]+=1}END{for(id in count)print id,count[id],term[id],class[id]}'|sort -k2nr,2 > go_count.out #4176

#合并go和term。（一个gene ID可能对应多列）
cat go_term.tab | awk 'BEGIN{FS=OFS="\t"}{plus=$3"("$4")";print $1,plus}' >go_plus_term.tab

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
END { if (rec) print rec }' > merge.awk

awk -f merge.awk go_plus_term.tab > go.anno.out
#得到最终注释结果

cp go.anno.out go_count.out ../result/













