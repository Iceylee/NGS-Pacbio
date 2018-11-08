#数字符个数
sed 's/[^"]//g' string | awk '{ print length }'


#替换文件名 破折号替换成下划线
for f in *.fasta
do 
mv -- "$f" "${f//-/_}"; 
done

#lists everything that isn't a directory (files, links, device files, etc.).
ls -p | grep -v / 



awk -f awkfile

#根据第2和5列，去重复
grep -v "\!" goa_uniprot_all.gaf |awk '!a[$2,$5]++' > goa_uniprot_all_duplicates_moved.gaf


#第三列
cat assoc.txt | cut -f 3 | head -20

#按照某列 匹配合并
#1join：两个文件 tab或者空格分隔都ok (70G 40min )
join -1 2 -2 2 -o 1.1 1.2 2.2 2.5 2.9 <(sort -k2,2 gene_swiss.id) <(sort -k2,2 goa_500.test) >result

#2awk：两个文件需要都是 tab分隔 (70G 10min)
#若2无重复，都是unique。或者需要输出unique的2列
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$5 FS $9;next}{print $0, a[$2]}' test.gaf gene_swiss.id >test.result
#若2列有重复。 小的表放前面，存入a。判断是否在大db中，有的话再输出
#将gene和swID存入a，然后db逐行判定。如果swID匹配则输出。
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$2]=$1;next}($2 in a){print a[$2],$0}' test.id test.gaf  >test.result




#blast结果转csv
cat cazy.bestHits|awk 'BEGIN{OFS=","}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'>cazy.csv

sed "1 i Gene_ID,cazy_ID,Identity,Align_Length,Mismatch,Gap,Q_start,Q_end,T_start,T_end,E_value,Score" cazy.csv

#sed
sed '1d' vfdbb.anno.out >> vfdba.anno.out #去首行 然后合并入文件a
cat vfdb-a.anno.out | sort -k1,1 -u>vfdb.anno.out  #按第一列排序并去重
#cat vfdba.anno.out vfdbb_no_head > vfdb.anno.out


cat 1b.txt | awk '{if (#~"00:")print($(x+1))}'|head

#pep,cds ID格式转为物种名+六位数字；locus起始和正负链
#>unitig_0|quiver_1 # 1 # 891 # 1 # ID=1_1;partial=10;start_type=Edge;rbs_motif=None;rbs_spacer=None;gc_cont=0.698
#>Pseudomonas_Putida_A316_000001	locus=1:891	+
BEGIN { FS = " # " ; OFS="\t" }
/>/ {	

	split($1,x,"|")
	split(x[2],y,"_")
	i = y[2]
	start = $2
	end = $3
	if ($4 == 1){
		chain = "+"
	}
	else{
		chain = "-"
	}
	printf (">Pseudomonas_Putida_A316_%06d\tlocus=%d:%d\t%s\n",i,start,end,chain)}
!/>/ {print $0}

#gff的第9列，简化只留物种名+六位ID
#跳过井号行
BEGIN {FS="\t";OFS="\t" }
!/\#/{
	split($9,x,";")
	split(x[1],y,"_")

	printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tID=pseudomonas_putidaA316_%08d\n",$1,$2,$3,$4,$5,$6,$7,$8,y[2])
}
/\#/{
	print($0)

}


#整列替换分隔符
cat cog2003-2014.csv |awk 'BEGIN{FS=",";OFS="\t"}{$1=$1}1'|head


#去除第一列末尾的分号
cat T3SS.txt |awk 'BEGIN{OFS="\t"}{if (NR!=1){print substr($1,1,length($1)-1)}}'|head

#打印9列以后所有列
awk '{ s = ""; for (i = 9; i <= NF; i++) s = s $i " "; print s }'

#计数
cat temp2 | awk 'BEGIN{OFS="\t"}{count[$15]+=1}END{for(id in count)print id,count[id]}'|sort -k2nr,2 >count.id

#将分号隔开的第二列分别输出。
#A0A1D6JJK6_MAIZE	Zm.17473;Zm.11111
#变成 A0A1D6JJK6_MAIZE	Zm.17473
#    A0A1D6JJK6_MAIZE Zm.11111
cat DownList_Protein.txt|awk 'BEGIN{OFS="\t"}{split($2,x,";");for (i=1;i<=length(x);i++)print $1,x[i] }' > temp1
#去除第二列为”“的row
cat temp1 |awk '{if($2!=""){print $0}}'>list
#取第二列 并去重复
cat Down.List|cut -f 2|sort -u > Down.Gene.List



#compare two 
comm -2 -3 <(sort file1) <(sort file2) > file3

#g:numerical  n:integer r:逆序
#第一列去重
sort -k1,1 -k12,12gr -k11,11g -k3,3gr $blastoutfile | sort -u -k1,1 --merge > bestHits 
#查看第一列是否sort了。完全sorted则不会返回任何结果
sort -k1,1 -c cog.bestHits


#打印匹配某patter（这里是多个大写字母）
#match的第三个参数y的[0]就是匹配到的字符串
cat cazy.bestHits | awk 'BEGIN{OFS="\t"}{split($2,x,"|");class=match(x[2],/([A-Z])+/,y);print y[0]}'|head

#条件判断 且
awk '$7=="+"&&FNR>1{OFS="\t";print $1,$4,$5}'

#处理多文件
cat rRNA.gff sRNA.gff tRNA.gff | awk 'FNR>1{print $1,$2,$3,$4}' >ncRNA.gff2
#FNR无效。cat已经将三个文件合并
awk 'FNR>1{print $1,$2,$3,$4}' rRNA.gff sRNA.gff tRNA.gff >ncRNA.gff2


#取起始4 + 最后10
cat MS_VS_BZ.sig.ID.blastp2kobas2identify.out| sed "1,4d"|tail -r|sed "1,10d"|tail -r > MSvsBZ.kegg.out




cat genome.fa.transdecoder.cds|awk '
BEGIN{FS=OFS=" "}
/>/{split($7,x,":");
s1=">" x[1];
s2=$4;
s3=$5;
gsub("[\\(\\)]","",$6);
s4="strand:" $6;
split(x[2],y,"(");
s5="pos:" y[1]
print s1,s2,s3,s4,s5
}
!/>/{print}' > novel.cds.fa

#大小写转换
tr '[:upper:]' '[:lower:]' < CNCI_nc.list > temp
cat test|awk '{split($1,x,"_");print tolower(x[1]) "_" x[2] "_" tolower(x[3])}'

#去掉大于号 -i即在原文件上修改
sed 's/>/ /g'
#
sed "1 i This is my monkey" pets.txt
#删除某行
sed "2d" pets.txt

#默认根据第一列来除？？
sort name2ID_from_gtf.list |uniq -u > name2ID_from_gtf_du.list

#匹配并打印 “XXX”
cat file|awk '{match($1,/\".*\"/,a);print a[0]}'|head
#
awk 'if not in' 等同与 in，不识别not，认为not是变量

#将一个文件按特定字符分成多个文件
awk '{print $0 "Sequence"> "file" NR}' RS='Sequence' out.dat

#将外部参数$name传给awk
awk -v var="$name" 'BEGIN{OFS="\t"}{print var,$0}'>$name.out

#grep -o match only
cat m6A.gff|grep "[0-9]+$" -o
cat m6A.gff|grep -E "[0-9]+$" -o #使用加号 需要E扩展grep

for i in $(find . -name "*.fastq");do gzip $i;done

for i in `ls *.fastq`
do
x=${i/.fastq/}
mv $i $x
done

for i in `ls ~/00.incipient_data/data_for_gene_prediction_and_RNA-seq/*.1.fastq`
do
i=${i/*\//}
i=${i/.1.fastq/}
done

#满足两种之一
cat q00000.keg|grep "^[AB]"|less
cat q00000.keg|grep -E  "^A[0-9]|^B  [0-9]"|less #使用|需要E扩展
cat temp.count |sed -r 's/^A[0-9]+//g'|head  #使用加号 用r扩展
#删除前3行
sed '1,3 d' file >fileout

#删除第3列
awk -i inplace '{$0=gensub(/\s*\S+/,"",3)}1' file


#awk match
match(string, regexp [, array])
 { if (match($0,/ your regexp/,m)) print m[0] }

cat ../assembly.fasta |seqkit grep -p unitig_18|seqret --filter -sbegin 1375 -send 1409

按string中的数字排序
sort -V 

#awk (not in array)
if (!(i in array))
###awk条件和语句要在同一行？ 留一个括号在条件行


#指定多个分隔符
cat AS.temp1|awk -F '[,/"]' 'BEGIN{OFS="\t"}{print $0,$2}' > AS.temp2

#替换分隔符 
cat ID.list|awk 'BEGIN{FS="\t";OFS=","}{$1=$1;print $0}' > temp
