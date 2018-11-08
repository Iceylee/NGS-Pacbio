# 0.1 gmap

# 除去注释，gmap比对上的数目，去重
grep -v -E "^@" gmap.sam|cut -f 1|sort -k1,1 -u|wc -l #245191
#flag 4为unmapped
grep -v -E "^@" gmap.sam|awk '$2==4{print $1}'|sort -k1,1 -u |wc -l #1513
#flag为0 or 16的reads 并去重
grep -v -E "^@" gmap.sam|awk '$2==0||$2==16{print $1}'|sort -k1,1 -u |wc -l #243678 #两者相加为总数245191
#认为大于10为唯一mapped
grep -v -E "^@" gmap.sam|sort -k1,1 -u|awk '$5>=10{n=n+1}END{print n}'   #220880 / 245191 = 90%


#质量分数XQ 为0 即比对到多个位置的,需要去重
grep -v -E "^@" gmap.sam|sort -k1,1 -u|awk '$5==0{n=n+1}END{print n}' #6935
#质量分数XQ 为40 即比对质量很好 $1==40 220311
#统计重复列个数
grep -v -E "^@" gmap.sam|awk '$5>=10{print $0}'|cut -f 1|sort |uniq -c
grep -v -E "^@" gmap.sam | awk -F"\t" 'BEGIN{print "flag\toccurrences"} {a[$2]++} END{for(i in a)print i"\t"a[i]}'


# 0.2 collapse
#总的一致性序列
seqkit stat isoform.fa #245191

#去冗余后的isoform条数
seqkit stat tofu.collapsed.rep.fa #70793
#平均每个PB collapse了几个reads
###计算字符个数 remained即参与collapsed的isoforms。
cat tofu.collapsed.group.txt|awk -F, '{sum = sum + NF}END {print sum}' #202874  / 70793 =2.86

#是否有重复 没有重复！！ 多个分隔符 ，换分隔符
cat tofu.collapsed.group.txt|awk -F'[\t,]' '
{$1=$1}
{for (i=2;i<=NF;i++){
print $i}}' > list
#gmap未比对上的 没有用
wc -l tofu.ignored_ids.txt #49701  245191-49701=


# 1.matchAnnot
##查看annotation.out最尾部
summary:   70791 isoforms read
summary:   70791 isoforms aligned, of which 0 were multiply mapped
summary:   69317 isoforms hit at least one gene, of which 757 were on opposite strand

##注释到已知基因的isoforms
69317


##不在基因区的isoforms
cat annotations.out|grep no_genes_found|wc -l
70791-69317=1474


Total isoform number 70791
#以下三项以参考序列为基准
Total unmap genes（种类）
Total mapped know genes
Total mapped know transcripts

#数据记录
#参考gtf的基因总数（种类）
cat my.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[2],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print a[0]}'|sort -k1,1 -u|wc -l
#49316
#


#isoform比对到参考gtf的基因种类。第三列original gene id
cat isoseq.info.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[3],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print a[0]}'|sort -k1,1 -u|wc -l
#18018-1=18017 （去掉no genes found）


#isoform未比对到参考gtf的。第二列 gene ID
cat ref.spec.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[2],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print a[0]}'|sort -k1,1 -u|wc -l #36684


#gtf有多少种tr ID  76486    第一列
awk 'BEGIN{FS=OFS="\t"}
($3=="transcript" || $3=="exon") {
split($9,x,";"); 
match(x[1],/\".*\"/,a);
gsub(/"/,"",a[0]);
print a[0]}' ../my2.gtf|sort -k1,1 -u|wc -l
#isoseq有多少种tr  24104
wc -l isoseq.tr.list
#ref.spec.gtf 52382   第一列
awk 'BEGIN{FS=OFS="\t"}
($3=="transcript" || $3=="exon") {
split($9,x,";");
match(x[1],/\".*\"/,a);
gsub(/"/,"",a[0]);
print a[0]}' ref.spec.gtf|sort -k1,1 -u|wc -l

#2.基因结构优化
wc -l structure_optimize.txt
 #39427-1
 #支持的reads数
cat structure_optimize.txt|cut -f 8|sort -k1,1 -u|wc -l
   #28396-1
cat extend_result.temp3|cut -f 3|sort -k1,1 -u|wc -l
   #16190


#3.lncRNA
wc -l angel.list #260
wc -l CNCI.list #1187
wc -l PLEK.list #877  
seqkit stat lncRNA.original.fasta #808
wc -l lncRNA.final.list #548 

seqkit stat novel.tr.fasta #926   (no_genes 1474=926+548)

#4.AS
#发生AS的基因个数
cat splice.as.out |cut -f 2|sort -k1,1 -u|wc -l #6124-1

#5.blast
wc -l *.anno.out
wc -l *.anno.out|grep -o -E "[0-9]+"|xargs -I {} echo {}/926|bc -l
wc -l pfam.anno.out #-1





#探讨gmap的MAPQ值
cat gmap.collapsed.sam|grep -E -o "XQ:.:[0-9]+"|grep -E -o "[0-9]+" > list

grep -v -E "^@" gmap.collapsed.sam|awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$17,$6}'|less -S


#gmap.collapsed.sam的第一列仍然有重复
#71063 （unqiue 70796）