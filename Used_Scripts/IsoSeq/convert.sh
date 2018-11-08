cat PB.list|awk '
BEGIN{FS=OFS="\t"}
{split($1,x,"|")
	print x[3],$1}' > read_PB.list

awk '
NR==FNR{a[$1]=$2;next}
/>/{split($1,x," ");split(x[1],y,">");p=y[2];print ">"a[p]" "x[2]}
!/>/{print $0}' ../../read_PB.list lncRNA.fasta > temp.fasta
mv temp.fasta lncRNA.fasta

awk '
NR==FNR{a[$1]=$2;next}
/>/{split($1,x," ");split(x[1],y,">");p=y[2];print ">"a[p]" "x[2]}
!/>/{print $0}' ../../read_PB.list novel.tr.fasta > temp.fasta
mv temp.fasta novel.tr.fasta


#blast
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2;next}
{print a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}
' ../../read_PB.list cog.anno.out > temp.anno.out
mv temp.anno.out cog.anno.out

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2;next}
{print a[$1],$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}
' ../../read_PB.list pfam.anno.out > temp.anno.out
mv temp.anno.out pfam.anno.out























#整理结果1: matchAnnot_result.txt
#1 of many中间都有空格，需要专门处理
cat annotations.out | grep "result"|awk '
BEGIN{FS=" ";OFS="\t"}
/1 of many/{print $2,$3" "$4" "$5" "$6,$7" "$8" "$9" "$10}
!/1 of many/{print $2,$3,$4}' > temp
# 去掉 " (1 of many)"
sed -- 's/ (1 of many)//g' temp > temp4matchannot
# 按照gtf 将name转换成id
#gene name - gene id - tr name - tr id
#use trID as index(to remove duplicates)
grep -v "#" ../Danio_rerio.GRCz10.86.gtf | awk '
BEGIN {FS=OFS="\t"}
{
	split($9,x,";")
	for (i in x) {
		if (x[i]~/gene_id/) {
			match(x[i],/\".*\"/,a)
			gsub(/"/,"",a[0])
			geneID=a[0];}
		if (x[i]~/transcript_id/) {
			match(x[i],/\".*\"/,b)
			gsub(/"/,"",b[0])
			trID=b[0];}
		if (x[i]~/gene_name/) {
			match(x[i],/\".*\"/,c)
			gsub(/"/,"",c[0])
			geneName=c[0];}
		if (x[i]~/transcript_name/) {
			match(x[i],/\".*\"/,d)
			gsub(/"/,"",d[0])
			trName=d[0];}
		}
	list[trID]=geneID"\t"geneName"\t"trID"\t"trName
}
END{
	for (i in list) print list[i]
}' >name2ID_from_gtf.list


#temp4matchanno按name2ID，将name转成ID
#temp5:geneLocus,geneName,trName,geneID,trID
awk '
BEGIN {FS=OFS="\t"}
NR==FNR {
	gene[$2]=$1
	tr[$4]=$3
	next
}
{
	if ($2=="no_genes_found")geneID = "no_genes_found"
	if ($2!="no_genes_found")geneID = gene[$2]
	trID = tr[$3]
	print $1,$2,$3,geneID,trID
}' name2ID_from_gtf.list temp > temp5
cat temp5|cut -f 1,4,5 | awk 'BEGIN{FS=OFS="\t"}$3!=""{print $0 ".mRNA"}$3==""{print}'|sed "1 i PacBio_transcript_id\tAnnoted_gene_id\tAnnoted_transcript_id"> matchAnnot_result.txt












