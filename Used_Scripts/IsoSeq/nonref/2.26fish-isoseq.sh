##############4.MatchAnnot(py2.7)########
# 0.1 gff to gtf
#cufflinks:gffread
#####gffread CID_femalev11_genemodel_edit.gff -T -o my.gtf
#??有少量序列报错

# 0.2 为避免matchanno报错，除去gtf中没有的染色体序列
ln -s ../ref/CID_annotation-ID.gtf my.gtf
cp ../gmap2/gmap.collasped.sorted.sam ./ #要修改 所以不用软连接
cat my.gtf |grep -v "#"|cut -f 1|sort -k1,1 -u > chr_in_gtf.list #32811
#sam文件的@行都去掉了，应该不影响？
awk 'NR==FNR{a[$1]=$1;next}($3 in a){print $0}' chr_in_gtf.list gmap.collasped.sorted.sam > my.sam 
#gtf文件有点问题：因此重新输出各列。
#直接match报错：ValueError: too many values to unpack
cat my.gtf|awk 'BEGIN{OFS="\t"}{print $1,$3,$4,$5,$6,$7,$8,$9,$11}' > my2.gtf


# 1.非冗余转录本比对结果与参考序列注释比较 
source activate Icey
matchAnnot.py --gtf my2.gtf my.sam > annotations.out
#遇到mismatching+长串字符
cat my.sam |grep "长串字符" #得到PB编号
cat my.sam |grep -v "PB.609.188|CI01000004:7620912-7632416(-)|i7_LQ_samplecb6195|c16763/f1p33/6070" > my2.sam

matchAnnot.py --gtf my2.gtf my2.sam > annotations.out



mkdir temp
cd temp
ln -s ../annotations.out ./



#整理结果1: matchAnnot_result.txt
#1 of many中间都有空格，需要专门处理
cat annotations.out | grep "result"|awk '
BEGIN{FS=" ";OFS="\t"}
/1 of many/{print $2,$3" "$4" "$5" "$6,$7" "$8" "$9" "$10}
!/1 of many/{print $2,$3,$4}' > temp

# 按照gtf 将name转换成id
#gene name - gene id - tr name - tr id
#use trID as index(to remove duplicates)
#修改，仅用transcript这一行
grep -v "#" ../my.gtf | awk '
BEGIN {FS=OFS="\t"}
{if($3 == "transcript"){
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
}}
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
cat temp5|cut -f 1,4,5 | sed "1 i PacBio_transcript_id\tAnnoted_gene_id\tAnnoted_transcript_id"> matchAnnot_result.txt
cp matchAnnot_result.txt ../result/


######
##很严格的gtf文件：包含gene name，gene id，tr name，tr id。且unique，name和id一一对应的话，直接用上面的方法。
##如果gtf是从gff文件转来，不满足上面的条件。可以直接写脚本修改gtf文件，使得gene name=gene id，tr name=tr id。
##再跑matchanno
##直接提取result内容。	

cat my2.gtf|awk '
BEGIN{FS=OFS="\t"}
{
	split($9,x,"; ")
	match(x[1],/\".*\"/,a)
	gsub(/"/,"",a[0])
	geneID = a[0]

	match(x[2],/\".*\"/,b)
	gsub(/"/,"",b[0])
	trID = b[0]

	col = x[1] "; " x[2] "; " "gene_name \"" geneID "\"; "    "transcript_name \"" trID "\"; " x[5]
	print $1,$2,$3,$4,$5,$6,$7,$8,col




}' |less -S

matchAnnot.py --gtf my3.gtf my2.sam > annotations.out

#接着上面的temp temp5 得到matchanno result



#整理结果2:isoseq.info.gtf 
#
cp ../../collapse/tofu.collapsed.gff ./

#trID PB.1.1 映射geneID
awk '
BEGIN{FS=OFS="\t"}
NR==FNR{split($1,x,"|");tr=x[1];ID[tr]=$2;next} 
{	split($9,x,";")
	match(x[2],/\".*\"/,a)
	gsub(/"/,"",a[0])
	qtr = a[0]
	if((qtr in ID)&& ID[qtr]!="no_genes_found") g=" orginal_gene_id \"" ID[qtr] "\"" ";"
	if((qtr in ID)&& ID[qtr]=="no_genes_found") g=" orginal_gene_id \"no_genes_found\";"
	if (qtr in ID ==0) g=""
	print $0 g

	}' matchAnnot_result.txt tofu.collapsed.gff > isoseq.info.gtf 
cp isoseq.info.gtf ../result/






#整理文件3：ref.spec.gtf
#原始gtf文件中所有转录本中没有被matchannot.result涵盖的转录本信息
#将isoseq有更正了的转录本信息，先从参考gtf中去除，之后换成isoseq的
#林霞：去除gene和CDS信息。transcript按ID，如果isoseq出现就去除该transcript。
#提取isoseq中有的transcript ID list
cat temp5|cut -f 5|sort -k1,1 -u > isoseq.tr.list
#处理原始gtf文件:第9行取第三个；匹配引号；去除引号；判断trID是否已在isoseq list中；无则输出

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {trList[$1]=1;next}
($3=="transcript" || $3=="exon") {
	split($9,x,";"); 
	match(x[3],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	if (!(a[0] in trList)) print $0}
' isoseq.tr.list ../my.gtf >ref.spec.gtf
cp ref.spec.gtf ../result/

##不去除gene和CDS信息
awk 'BEGIN{FS=OFS="\t"}
NR==FNR {trList[$1]=1;next}
{
	split($9,x,";"); 
	match(x[3],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	if (!(a[0] in trList)) print $0}
' isoseq.tr.list ../my.gtf >ref.spec.gtf2
cp ref.spec.gtf ../result/



#整理文件4：merge.gtf
#isoseq.info.gtf + ref.spec.gtf
cat ref.spec.gtf isoseq.info.gtf > merge.gtf
cp merge.gtf ../result/


