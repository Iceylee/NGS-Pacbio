##############1.combine all hq&lq,and get fasta file##############
mv *hq* hq.fq
mv *lq* lq.fq
seqret -sequence hq.fq -outseq hq.fa
seqret -sequence lq.fq -outseq lq.fa
cat hq.fa lq.fa > isoform.fa
seqkit stat *.fa



ref_fa=GCF_000005005.2_B73_RefGen_v4_genomic.fna
ref_name=Maize_B73

##############2.gmap比对参考基因组########
#1.参考基因组索引  #1.3G 40min  #1h 900M
gmap_build -d $ref_name ref/$ref_fa

#2.比对  -D 索引所在目录 -d db前缀同文件夹名字 
#-t 20 500M 30min
#-d gmap_db
mkdir gmap && cd gmap
nohup gmap -D /home/liyb/miniconda3/share/$ref_name -d $ref_name -f samse -n 0 -t 20 ../raw/isoform.fa > gmap.sam 2> gmap.sam.log &

#3.sort
sort -k 3,3 -k 4,4n gmap.sam > gmap.sorted.sam

##############3.ToFu########
source activate anapy2.7.9
export VENV_TOFU=~/biosoft/VENV_TOFU
source $VENV_TOFU/bin/activate

#1.融合基因
cd Fusion
fusion_finder.py --input ../raw/isoform.fa -s ../gmap/gmap.sorted.sam -o fusion 
#--cluster_report_csv cluster_report.csv cDNA-primer版本没有这个参数。csv在SMRTLink job 里面得到。

#2.非冗余转录本
cd collapse
collapse_isoforms_by_sam.py --input ../raw/isoform.fa -s ../gmap/gmap.sorted.sam --dun-merge-5-shorter -o tofu
deactivate
source deactivate 

#3.去冗余后重新比对和sort 
cd ../gmap2
nohup gmap -D /home/liyb/miniconda3/share/$ref_name -d $ref_name -f samse -n 0 -t 20 ../collapse/tofu.collapsed.rep.fa > gmap.collapsed.sam 2> gmap.collapsed.sam.log &
sort -k 3,3 -k 4,4n gmap.collapsed.sam > gmap.collasped.sorted.sam



#整理文件
mv gmap.collasped.sorted.sam collapsed.sam
gzip collapsed.sam -c > collapsed.sam.gz
gzip gmap.sorted.sam -c > gmap.sorted.sam.gz



##############4.MatchAnnot(py2.7)########
# 0.1 gff to gtf
#cufflinks:gffread
gffread *.gff -T -o my.gtf
#??有少量序列报错

# 0.2 为避免matchanno报错，除去gtf中没有的染色体序列
ln -s ../ref/my.gtf my.gtf
ln -s ../gmap2/gmap.collasped.sorted.sam ./
cat my.gtf |grep -v "#"|cut -f 1|sort -k1,1 -u > chr_in_gtf.list 
#sam文件的@行都去掉了，应该不影响？
awk 'NR==FNR{a[$1]=$1;next}($3 in a){print $0}' chr_in_gtf.list gmap.collasped.sorted.sam > chr_filter.sam 
#!检查sam文件的行数 确保相差不大

#gtf文件
#报错：tranName = re.search(regexTran, attrs).group(1)
#AttributeError: 'NoneType' object has no attribute 'group'
#表示需要transcriptName

cat my.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[1],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print $0 " transcript_name \""a[0] "\";"} ' > my2.gtf

#exon number



# 1.非冗余转录本比对结果与参考序列注释比较 
source activate Icey

#加参数--format alt 不报错“缺少exon信息”
matchAnnot.py --gtf my2.gtf --format alt chr_filter.sam > annotations.out




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
grep -v "#" ../my2.gtf | awk '
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
	list[trID]=geneID"\t"trID"\t"geneName"\t"trName
}
END{
	for (i in list) print list[i]
}' >name2ID_from_gtf.list

#调整gtf第9列顺序：geneID trID geneName trName
cat my2.gtf|awk '
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
	$9 = "gene_id " "\"" geneID "\""     "; transcript_id "  "\"" trID "\""          "; gene_name "  "\"" geneName "\""      "; transcript_name "    "\"" trName "\""  ";"
	print $1,$2,$3,$4,$5,$6,$7,$8,$9}' >my3.gtf

#temp4matchanno按name2ID，将name转成ID
#temp5:geneLocus,geneName,trName,geneID,trID

#特殊：这个gtf文件里，同一个gene name对应不同的gene ID。而tr name其实就是tr ID。因此全部以tr ID 来对应

awk '
BEGIN {FS=OFS="\t"}
NR==FNR {
	gene[$2]=$1
	tr[$2]=$4
	next
}
{
	if ($2=="no_genes_found")geneID = "no_genes_found"
	if ($2!="no_genes_found")geneID = gene[$3]
	trID = tr[$3]
	print $1,geneID,trID,$2,$3
}' name2ID_from_gtf.list temp > temp5
cat temp5|cut -f 1,4,5 | sed "1 i PacBio_transcript_id\tAnnoted_gene_id\tAnnoted_transcript_id"> matchAnnot_result.txt




#整理结果2:isoseq.info.gtf 
#补充origninal结果
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
 
#check??应该可以一一对应，但因为matchanno自身的bug，输入的是过滤了少量染色体的sam文件。因此会有稍微的差异
cat isoseq.info.gtf awk '$3=="transcript" {print $0}' |wc -l 
   #58097 多于matchAnnot，因为少量chr没有在matchannot中#not found：26
wc -l matchAnnot_result.txt 
   #58073-1（首行）
cat temp9|grep "not found"|awk '$3=="transcript" {print $0}'|wc -l #26
cat temp9|grep "orginal_gene_id"|awk '$3=="transcript" {print $0}'|wc -l #57231
cat temp9|grep "no_genes"|awk '$3=="transcript" {print $0}'|wc -l #840







#整理文件3：ref.spec.gtf
#原始gtf文件中所有转录本中没有被matchannot.result涵盖的转录本信息
#将isoseq有更正了的转录本信息，先从参考gtf中去除，之后换成isoseq的
#林霞：去除gene和CDS信息。transcript按ID，如果isoseq出现就去除该transcript。
#提取isoseq中有的transcript ID list
cat temp5|cut -f 5|sort -k1,1 -u > isoseq.tr.list
#处理原始gtf文件:第9行取第二个trID（确认一下）；匹配引号；去除引号；判断trID是否已在isoseq list中；无则输出

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {trList[$1]=1;next}
($3=="transcript" || $3=="exon") {
	split($9,x,";"); 
	match(x[2],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	if (!(a[0] in trList)) print $0}
' isoseq.tr.list ../my3.gtf >ref.spec.gtf2








#整理文件4：merge.gtf
#isoseq.info.gtf + ref.spec.gtf
cat ref.spec.gtf2 isoseq.info.gtf > merge2.gtf


#result
mkdir result
cp ../temp/merge.gtf ./
cp ../temp/ref.spec.gtf ./ 
cp ../temp/isoseq.info.gtf ./
cp ../temp/matchAnnot_result.txt ./



#check
#gtf特有的tr ID 种类数
cat ref.spec.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[3],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print a[0]} ' |sort -k1,1 -u >ref.spec.tr.list
cat ref.spec.tr.list|wc -l #45787
#isoseq有的tr ID种类数
wc -l isoseq.tr.list #12488
#求交集
comm ref.spec.tr.list isoseq.tr.list >temp.comm
cat temp.comm|awk '{print $3}'|sort -k1,1 -u|wc -l #1 仅有一个空值，表明没有交集








#############lncRNA
ln -s ../MatchAnno/result/matchAnnot_result.txt
ln -s ../collapse/tofu.collapsed.rep.fa
cat matchAnnot_result.txt |grep "no_genes_found"|cut -f 1 > no_genes.list
#提取相应序列
cat tofu.collapsed.rep.fa |seqkit grep -f no_genes.list >nc.fa

#PELK (py2.7)
source activate Icey
PLEKDIR=/home/liyb/biosoft/PLEK.1.2
python $PLEKDIR/PLEK.py -fasta nc.fa -out PLEK.out -thread 20 1>PLEK.log
cat PLEK.out |grep Non-coding >PLEK_nc.list  

#CNCI(py2.7) 慢 5h

CNCIDIR=/home/liyb/biosoft/CNCI-master 
#-m 脊椎动物或者植物 参数 ve pl
nohup python $CNCIDIR/CNCI.py -f nc.fa -o CNCI_out -p 15 &
cat CNCI_out/CNCI.index |awk 'BEGIN{OFS="\t"}$2=="noncoding"{print $1,$2}' > CNCI_nc.list

#合并
mv *nc.list result/
#提取gene list
#不同文库的情况

cat PLEK_nc.list |awk '{split($3,x," ");print x[1]}'|sed 's/>//g' >PLEK.list
cat CNCI_nc.list |awk 'BEGIN{FS=" "}{print $1}'> CNCI.list

#取交集
comm -12 <( sort CNCI.list ) <( sort PLEK.list ) >CNCI_PLEK.list 

##去掉比对上pfam库的序列
ln -s ../../blast/result/pfam.anno.out
cat pfam.anno.out|cut -f 1 > pfam.list
#仅在CNCI_PLEK_CPAT.list 中出现的序列
comm -23 <( sort CNCI_PLEK.list ) <( sort pfam.list ) >CNCI_PLEK_pfam.list 

#提取lncRNA序列
cat ../tofu.collapsed.rep.fa |seqkit grep -f CNCI_PLEK_pfam.list  > lncRNA.original.fasta

#ORF预测 有编码潜能的序列
mkdir angel && cd angel
source activate Icey
dumb_predict.py ../lncRNA.original.fasta angel.dumb --cpus 15
#--min_aa_length 300 

cd ..
cat angel.dumb.final.cds|grep ">" |cut -d " " -f 1|sed 's/>//g'|sed 's/|m.*//g' > angel.list
#去掉angel.list中的序列(23条)
cat lncRNA.original.fasta |seqkit grep -f angel.list -v > lncRNA.fasta


#新转录本 262-23=239
cat lncRNA.fasta|grep ">" |cut -d " " -f 1 |sed 's/>//g' > lncRNA.final.list
#取差集 840-239=601
comm -23 <( sort ../no_genes.list ) <( sort lncRNA.final.list ) >novel.tr.list 
#提取序列
cat ../tofu.collapsed.rep.fa |seqkit grep -f novel.tr.list > novel.tr.fasta



#############blastx
#nr,trembl,pfam
#nr
cd ~/analysis/pro14_Maize/blast
cDNA=novel.tr.fasta
db=nr #change

dbdir=${db}dir
dbdir=${!dbdir}
log=$db.log
out=$db.out
echo $dbdir $cDNA $log $out

nohup blastx -db $dbdir -query $cDNA -num_threads 30 -out $out -outfmt "6" -evalue 1e-5 -max_target_seqs 1 > $log 2>&1 &


