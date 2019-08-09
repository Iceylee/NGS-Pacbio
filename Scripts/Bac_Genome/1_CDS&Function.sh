###去掉|quiver
cat assembly.fasta|awk '/>/{split($1,x,"|");print x[1]}!/>/{print $0}' > temp.fasta
mv temp.fasta assembly.fasta

mkdir prodigal
mkdir blast
mkdir blast/cog
mkdir blast/pfam
mkdir blast/GO
mkdir blast/bestHits
mkdir blast/result
mkdir blast/out
mkdir blast/kegg

#############0.调整ORF（若单contig）#############
#name should be assembly.fasta
python ~/biosoft/oriFinder/oriFinder.py assembly.fasta

#调整起始点
#255794
cat assembly.fasta |seqret --filter -sbegin 1 -send 255794 > 1.fasta
cat assembly.fasta |seqret --filter -sbegin 255794 -send 4467693 > 2.fasta
cat 2.fasta 1.fasta |union -filter > new.fasta

#############1.CDS (prodigal)#############
cd prodigal
name=Brevibacillus_sp.WF146

prodigal -i ../assembly.fasta -d prodigal.cds -f gff -o prodigal.gff -s prodigal.stat -a prodigal.pep
#prodigal -i ../assembly.fasta -o prodigal.gbk
#-f指定输出gff格式  默认为gbk

#1contig选择gff1和pep1；多个contig选择gff2和pep2
#pep输出名字修改:物种名字+六位编号  locus和正负链
scriptdir=~/analysis/scripts
pepawk=pep1.awk
gffawk=gff1.awk

#pep
awk -v name=$name -f $scriptdir/$pepawk prodigal.pep >${name}.pep
#cds
awk -v name=$name -f $scriptdir/$pepawk prodigal.cds >${name}.cds
#gff 提取最后一列的ID，并修改名字。
awk -v name=$name -f $scriptdir/$gffawk prodigal.gff >${name}.gff

cp ${name}.pep ../assembly.faa
################2.blast#############
#共10个文件，按顺序
#nr,trembl,pfam
#sp,cog


cd /home/liyb/analysis/pro8_0001_1contig/blast
pep=/home/liyb/analysis/pro8_0001_1contig/assembly.faa
db=cog #change

dbdir=${db}dir
dbdir=${!dbdir}
log=$db.log
out=$db.out
echo $dbdir $pep $log $out

nohup blastp -db $dbdir -query $pep -num_threads 20 -out $out -outfmt "6" -evalue 1e-5 -max_target_seqs 1 > $log 2>&1 &

#nr分成两部分跑的
cat nr.half.out nr.1.out >> nr.out
sort -k1,1 -u nr.out > ../bestHits/nr.bestHits

ls -sh
ls -alt

################2.get best hits
#blast用了(-maxtargetseqs参数)
##Gene_ID\tSwissProt_ID\tIdentity\tAlign_Length\tMismatch\tGap\tQ_start\tQ_end\tT_start\tT_end\tE_value\tScore
#based on score
#mv {}.out out/ #blast完毕转入out文件夹
 
for i in "sp" "nr" "cog" "trembl" "pfam"
do
	sort -k1,1 -u out/${i}.out > bestHits/${i}.bestHits
done 



#blast未用参数-maxtargetseqs
#sort -k1,1 -k12,12gr -k11,11g -k3,3gr $blastoutfile | sort -u -k1,1 > bestHits

################3.anno
#根据geneid在anno中提取注释信息，作为最后一列加入，得到添加注释列的文件.anno.out
#空格分隔，读取后，输出第一列，后面所有列的合并作为第二列。
annodir=~/analysis/database/anno
for i in "sp" "nr" "cog" "trembl"
do
	python ~/analysis/scripts/db_anno_extract.py $annodir/${i}.anno bestHits/${i}.bestHits result/${i}.anno.out
done 

mv result/cog.anno.out result/cog.anno.out.orignal

################4.需要特殊处理才能得到的.anno.out
#1.pfam
cd pfam
ln -s ../bestHits/pfam.bestHits
ln -s ~/analysis/database/Pfam-A.regions.uniprot.tsv ./
ln -s ~/analysis/database/pfam_term_dup_moved.tab ./

cat pfam.bestHits |awk 'BEGIN{OFS="\t"}{split($2,x,"_");print $1,x[1]}'>gene-swiss.id
#gene-swissprot-pfamID
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$5;next}{print $0, a[$2]}' Pfam-A.regions.uniprot.tsv gene-swiss.id >temp
#pfamID匹配term
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{print $0, a[$3]}' pfam_term_dup_moved.tab temp > pfam.out
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$3;b[$1]=$4;next}{print $1, a[$1],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,b[$1]}' pfam.out pfam.bestHits > pfam.anno.out
cp pfam.anno.out ../result/

#2.cog
cd cog
ln -s ~/database/cog/fun2003-2014.tab ./
ln -s ~/database/cog/cognames2003-2014.tab ./
ln -s ~/analysis/database/cog2003-2014.csv.tab ./
ln -s ../result/cog.anno.out.orignal ./
#key第一列，需要第7列； 按照第二列的 gi|91777738|ref|YP_552946.1| 第二个查询
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$7;next}{split($2,x,"|");y=x[2];print $0, a[y]}' cog2003-2014.csv.tab cog.anno.out.orignal >temp
#key第一列，需要第二列和第三列；按照14列COG编号查询
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2 FS $3;next}{print $0, a[$14]}' cognames2003-2014.tab temp > temp2



#按照第15列分组计数。用于作图 (第二列，按数字，从大到小倒序排列)
cat temp2 | awk 'BEGIN{OFS="\t"}{count[$15]+=1}END{for(id in count)print id,count[id]}'|sort -k2nr,2 >count.id
#多个类的拆分开，并加入各类的计数
cat count.id | awk '
	BEGIN{OFS="\t"}
	{
		split($1,x,"")
		for (i in x){
			class = x[i]
			count[class]+=$2
		}
	}
	END{for (id in count)print id,count[id]}' |sort -k1,1 > count.split.id





#count id加上一列，说明分类的意义；key第一列，需要第二列；按照第一列查询
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{print $0, a[$1]}' fun2003-2014.tab count.split.id > cog_count_fun.id
#合并最后几列 分隔符 ”;“
cat temp2 | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13";"$14";"$15";"$16}'> temp3
#cog coding conversion
iconv -f Windows-1252 -t UTF-8//TRANSLIT temp3 > temp4
mv temp4 cog.anno.out

cp cog.anno.out cog_count_fun.id ../result/



#3.kegg
##下载q00000.keg，query.ko.txt
#有kegg注释的 
cat query.ko.txt |awk '$2!=""{i=i+1}END{print i}' #2508 每个基因unique
less -S q00000.keg |grep "^D"|grep -E ".+; K" -o|sort -k2,2 -u|wc -l #2508 每个基因可能有多个KO对应，因此去重复
cat q00000.keg|grep "^D"|wc -l #4878 不去重复 
#提取D行，去掉第一列（D空格空格）。按分号+空格分隔，并将kegg与后面term之间的空格转成冒号
cat q00000.keg |grep "^D"|sed 's/^D      //g'|awk 'BEGIN{FS="; ";OFS="\t"}{sub(" ",":",$2);print $1,$2";"$3}' |sort -k1,1 >temp.anno #4878
#很多重复（去重复）
cat temp.anno|awk '
BEGIN{OFS="\t"}{split($2,x,":");print $1,$2,x[1]}'|awk 'BEGIN{OFS="\t"}!a[$1,$3]++'|cut -f 1,2 >kegg.anno.out


#取出所有信息，按A，B D 格式
cat q00000.keg |awk 'BEGIN{OFS="\t"}
/^A/{A_name=$0}
/^B  /{B_name=$0}
/^D/{D_name=$0;print A_name,B_name,D_name}' |awk '{OFS="\t"}{split($3,x,";");split(x[1],y,"_");print $1,$2,$3,y[3]}' > temp1
#去除KO，B分类完全相同的行（为何3排序可以去除？还不清楚）
cat temp1|sort -k3 -u > temp_du

cat temp_du|awk 'BEGIN{OFS="\t"}{a[$2]++;b[$2]=$1}END{for (i in a)print b[i],i,a[i]}' > temp.count
cat temp.count|sed -r 's/A[0-9]+ //g'|sed -r 's/B  [0-9]+ //g'| sort -k3,3 -n >kegg.count #sum 
cp kegg.count kegg.anno.out ../result/

##提取pathway ID
awk 'BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$3"\t"$4;next}
{split($2,x,":");
K_ID=x[1];
term=a[K_ID];
print $1,$2,K_ID,term}' ko00000_kegg.2018.5.30.txt kegg.anno.out >temp

cat temp|awk 'BEGIN{FS=OFS="\t"}
{match($4,/ko[0-9]+/,a);print $1,$2,$3,a[0],$5}' > kegg_pathway.out



##

################5.其他需要的文件
#1.genome的gene_list

cat ../assembly.faa |grep ">" > temp.list
sed -i 's/>//g' temp.list #在原文件删除 匹配大于号删除
cat temp.list|awk 'BEGIN{FS=" "}{print $1}' > result/genome.list #去除后面的locus内容
#2.每个contig长度和GC含量
#GC content
seqkit fx2tab -l -g -n -i -H assembly.fasta > test.tab
#3.转格式
###blast_info.R 输出out文件夹
###
cd result/
python ~/work/scripts/txt2xls.py out/
###
python ~/work/scripts/txt2xls.py txt/








