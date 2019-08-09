#############基因结构优化
#提取iso，gene，tr信息，合并到exon栏 
mkdir optimization
cd optimization

ln -s ../MatchAnno/annotations.out
cat annotations.out |awk '
BEGIN{OFS="\t"}
$0~/isoform:/ {split($0,x," ");iso=x[2]}
$0~/gene:/ {split($0,x," ");gene=x[2]}
$0~/tr:/{split($0,x," ");tr=x[2]}
$0~/exon:/ {split($0,x," ");print iso,gene,tr,x[4],x[5],x[6],x[7],x[8],x[9]}' > iso_gene_tr_exon.tab
# 去掉起始 终止带 句号；即没有匹配上的
cat iso_gene_tr_exon.tab|awk '($4!="."&&$5!="."&&$7!="."&&$8!="."){print $0}' >exon.temp1
# 去掉均为0即完全匹配的
cat exon.temp1|awk '($6!="0"||$9!="0"){print $0}' > exon.temp2
# 每个转录本，只保留最优result对应的tr ID;并对应gene ID	
ln -s ../MatchAnno/temp/temp5
awk  'BEGIN{OFS="\t"}
	NR==FNR{a[$1]=$3;b[$2]=$4;next}
	a[$1]==$3{print $0,b[$2]}
' temp5 exon.temp2 > exon.temp3 
#只保留两端至少一端可以延长的
cat exon.temp3|awk '($6>0)||($9<0){print $0}' > exon.temp4


#整理
#从$1提取strand和chr
cat exon.temp4|awk '
BEGIN{OFS="\t"}
	{split($1,x,"|");
	split(x[2],y,":");
	chr=y[1];
	match(x[2],/\(.\)/,a);
	gsub(/[\(\)]/,"",a[0]);
	strand=a[0];
	print $0,chr,strand}' > extend_result.temp1

# 坐标整理：original(tr);extended;final(iso)
# cat extend_result.temp1|awk '
# 	BEGIN{OFS="\t"}
# 	$6>0&&$9<0{strand="3,5";print $0,$5"-"$8,$4"-"$7,$4"-"$7,strand}
# 	$6>0&&$9>=0{if($12=="+")strand="5";if($12=="-")strand="3";print $0,$5"-"$8,$4"-"$7,$4"-"$8,strand}
# 	$6<=0&&$9<0{if($12=="+")strand="3";if($12=="-")strand="5";print $0,$5"-"$8,$4"-"$7,$5"-"$7,strand}' > extend_result.temp2
# #输出列名，调整顺序
# cat extend_result.temp2|awk '
# BEGIN{OFS="\t";print "Gene_ID\t3_or_5_end\tchr\tStrand\tOriginal_region\tIsoform_region\tFinal_region\tReads_id"}
# {print $10,$16,$11,$12,$13,$14,$15,$1}' > extend_result.out
# mv extend_result.out structure_optimize.txt


#报告整理II
###将isoform坐标改为延长区域坐标。3，5均延长的分成两行表示
# 原始顺序：4iso起始 5tr起始 6差异 7iso终止 8tr终止 9差异
cat extend_result.temp1|awk '
	BEGIN{OFS="\t"}
	$6>0&&$9<0{
		if($12=="+"){
			strand="3";
			print $0,$5"-"$8,$8"-"$7,$4"-"$7,strand
			NR=NR+1
			strand="5"
			print $0,$5"-"$8,$4"-"$5,$4"-"$7,strand
		}
		if($12=="-"){
			strand="3";
			print $0,$5"-"$8,$4"-"$5,$4"-"$7,strand
			NR=NR+1
			strand="5"
			print $0,$5"-"$8,$8"-"$7,$4"-"$7,strand
		}

	}
	$6>0&&$9>=0{
		if($12=="+")strand="5";
		if($12=="-")strand="3";
		print $0,$5"-"$8,$4"-"$5,$4"-"$8,strand
		}
	$6<=0&&$9<0{
		if($12=="+")strand="3";
		if($12=="-")strand="5";
		print $0,$5"-"$8,$8"-"$7,$5"-"$7,strand
		}' >extend_result.temp3

cat extend_result.temp3|awk '
BEGIN{OFS="\t";print "Gene_ID\t3_or_5_end\tChr\tStrand\tOriginal_region\tExtended_region\tFinal_region\tReads_id"}
{print $10,$16,$11,$12,$13,$14,$15,$1}' > extend_result.out
mv extend_result.out structure_optimize.txt
#每行表示一个修正的外显子。因此同一个gene+isoform，可以有多行。因为有多个外显子。


#############2.可变剪切 
mkdir AS
cd AS

ln -s ../MatchAnno/result/merge.gtf
# 需要transcript id按同样的顺序。比如都是9列第二个
# 原始gtf的9列第二个是gene_version.去除并只输出geneid和trid
# cat merge.gtf|awk '
# BEGIN{FS=OFS="\t"}
# $9~/gene_version/{
# 	split($9,x,";");
# 	$9=x[1]";"x[3]";"
# 	print $0
# }
# $9!~/gene_version/{
# 	print $0
# }' >temp.gtf


#输出9列的第1个和第2个
cat merge.gtf|awk '
BEGIN{FS=OFS="\t"}
{
split($9,x,";");
$9=x[1]";"x[2]";"
print $0
}' >temp.gtf

#报错的行：第二列 Curated Genomic 有空格。手动合并成CuratedGenomic


# Pairwise AS event extraction from a transcriptome 成对
astalavista -t asta -i temp.gtf #as1.gtf 105305
# Complete AS events
# 如果有三个逗号的，这里会作为一列显示。而前者会分成两组显示
#astalavista -t asta -i temp.gtf -d 0 #as2.gtf 49002


mv temp_sorted.gtf_astalavista.gtf as.gtf
cat as.gtf|awk '
BEGIN{OFS="\t"}
$9~/structure \"0,1-2\^\"/ {
	print $0,"Exon_Skipping"
}
$9~/structure \"1\^,2\^\"/ {
	print $0,"Alternative_Donors"
}
$9~/structure \"1-,2-\"/ {
	print $0,"Alternative_Acceptors"
}
$9~/structure \"0,1\^2-\"/ {
	print $0,"Intron_Retention"
}
$9~/structure \"1-2\^,3-4\^\"/ {
	print $0,"Mutually_Exclusive_Exons"
}
$9!~/structure \"0,1-2\^\"/ && $9!~/structure \"1\^,2\^\"/ && $9!~/structure \"1-,2-\"/ && $9!~/structure \"0,1\^2-\"/ && $9!~/structure \"1-2\^,3-4\^\"/ {
	print $0,"Other"
}
'  > AS.temp1

# tr ID匹配gene ID
# AS的结果，用逗号分割再用/，随便取一个就行。取第一个tr ID，添加为最后一列(11)
cat AS.temp1|awk -F '[,/"]' 'BEGIN{OFS="\t"}{print $0,$2}' > AS.temp2
# 得到tr ID 对应gene ID tab（部分来自gtf，部分来自temp5）
ln -s ../MatchAnno/temp/name2ID_from_gtf.list
ln -s ../MatchAnno/temp/temp5
# 来自原始gtf的geneID对应trID
cat name2ID_from_gtf.list |cut -f 1,3 > gene_tr.gtf
# 来自matchanno的geneID对应isoID（PB1.1）
cat temp5|awk '
BEGIN{OFS="\t"}
{split($1,x,"|");print $4,x[1]}' > gene_tr.iso
# 合并
cat gene_tr.gtf gene_tr.iso > gene_tr.tab
#匹配trID
awk '
BEGIN{OFS="\t"}
NR==FNR{a[$2]=$1;next}
{print $0,a[$11]}
	' gene_tr.tab AS.temp2 > AS.temp3

# PB1.1 匹配 PB1.1|1：111-222（+）|c111/f1p1/111
ln -s ../collapse/tofu.collapsed.rep.fa
cat tofu.collapsed.rep.fa |grep ">"|sed 's/>//g' > PB.list
cat PB.list|awk '
	BEGIN{OFS="\t"}
	{split($1,x,"|");print x[1],$1}' > PB_total.list

#提取结果
cat AS.temp3|awk '
BEGIN{FS=OFS="\t"}
{	split($9,x,"; ")
	match(x[1],/\".*\"/,a)
	gsub(/"/,"",a[0])
	split(a[0],x,",")
	match(x[4],/\".*\"/,b)
	gsub(/"/,"",b[0])
	print $1,$12,$10,$4,$5,x[1],x[2]}' > AS.temp4

#PB替换成更长的表达，用分号空格作为分隔。如果只有一个，去掉分隔符。
cat AS.temp4|awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$2;next}
{	PBlist=""
	split($6,x,"/")
	for (i in x){
		if(match(x[i],/PB\.[0-9]+\.[0-9]+/,b)){
			oldPB=b[0]
			newPB=a[oldPB]
			PBlist=newPB"; "PBlist		
		}
		else {
			PBlist=x[i]"; "PBlist
		}
	}

	PBlist2=""
	split($7,x,"/")
	for (i in x){
		if(match(x[i],/PB\.[0-9]+\.[0-9]+/,b)){
			oldPB=b[0]
			newPB=a[oldPB]
			PBlist2=newPB"; "PBlist2		
		}
		else {
			PBlist2=x[i]"; "PBlist2
		}
	}
	print $1,$2,$3,$4,$5,PBlist,PBlist2
	
}' PB_total.list AS.temp4  >AS.temp5

cat AS.temp5|sed "1 i Chromosome\tGene_ID\tAS_type\tEvent_start\tEvent_end\tTranscript_ID1\tTranscript_ID2" > splice.as.out


# count
cat splice.as.out|awk '{i[$3]++} END{for (event in i)print event,i[event]}' >event_count.out
cat splice.as.out|awk '
$3=="Other"{a[$2]=$2}
$3=="Alternative_Donors"{b[$2]=$2}
$3=="Alternative_Acceptors"{c[$2]=$2}
$3=="Intron_Retention"{d[$2]=$2}
$3=="Exon_Skipping"{e[$2]=$2}
$3=="Mutually_Exclusive_Exons"{f[$2]=$2}
END{print "Other: "length(a)"\nAlternative_Donors: "length(b)"\nAlternative_Acceptors: "length(c)"\nIntron_Retention: "length(d)"\nExon_Skipping: "length(e)"\nMutually_Exclusive_Exons: "length(f)}'>gene_count.out

#check
cat splice.as.out|grep "Mutually"|sort -k2,2 -u|wc -l




#############饱和曲线
mkdir Curve && cd Curve


cp /home/software/prog/smrtlink_data/jobs/000/000388/tasks/pbtranscript.tasks.separate_flnc-0/combined/all.cluster_report.csv ~/analysis/pro14_Maize/Curve/report.csv

# 多个CCS对应到1个cluster（isoform）
sed '1 d' report.csv |grep -v "NonFL"|awk 'BEGIN{FS=",";OFS="\t"}{print $2,$1}' > read_cluster.tab
# 多个isoform（tr）对应到1个gene
ln -s ../MatchAnno/result/matchAnnot_result.txt
cat matchAnnot_result.txt |grep -v "no_genes_found" > matchAnnot_result.genes




#smrt analysis，sequel下机数据
# 整合（cluster即isoform）



#cluster ：i1_ICE_sample4f5084|c15221  gene:i1_HQ/LQ_sample4f5084|c15221
cat read_cluster.tab |awk '{split($2,x,"|");print substr(x[2],2)}'|sort -k1n,1 -u > cluster.list #awk去除字符c 从第二位开始读取

#验证cluster是独一无二的，无论i的编号
##得到i2_LQ_sample4f5084|c11716 （后面的删掉）
cat matchAnnot_result.genes|awk -F '[/|]' '{print $3"|"$4}'|sort -k1,1 -u |wc -l  #awk 指定多个分隔符
#结论：每个i都从c0开始编号，因此不是unique的

#验证同一个i：cluster，不会即是hq又是lq
##删除HQ和LQ，然后去重计数
cat matchAnnot_result.genes|awk -F '[/|]' '{print $3"|"$4}' > quest1
cat quest1|awk -F '_' '{print $1"_"$3}'|sort -k1,1 -u|wc -l
wc -l quest1
#结论：hq与lq可删除。每一个cluster要么对应hq，要么对应lq


###########

#1.matchAnnot_result.genes获取i2_LQ_sample4f5084|c11716，并去掉LQ与HQ
cat matchAnnot_result.genes|awk -F '[/|\t]' '{print $3"|"$4"\t"$7"\t"$8}' > temp1
cat temp1|awk -F '_' '{print $1"_"$3}'>temp2

#2.去除read_cluster.tab中第二列的ICE
cat read_cluster.tab|awk -F '\t' '{split($2,x,"_");print $1"\t"x[1]"_"x[3]}'>temp3

mv temp2 isoform_gene.tab
mv temp3 read_cluster.tab


#smrt portal
sed '1 d' matchAnnot_result.genes|awk '
BEGIN{FS=OFS="\t"}
{split($1,x,"|");split(x[3],y,"/");print y[1],$2}' > isoform_gene.tab


#portal & analysis
awk '
BEGIN{FS="\t";OFS="\t"}
NR==FNR{a[$1]=$2;next}
{print $1,$2,a[$2]}
' isoform_gene.tab read_cluster.tab > read_isoform_gene.tab

#5000为单位 计数
cat read_isoform_gene.tab|awk '
BEGIN{OFS="\t";print"FL_reads_number\tIsoform_number\tGene_number"}
{i++;a[$2]=$2;b[$3]=$3;if (i%5000==0)print i,length(a),length(b)}' > saturation_curve.txt








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

#CPAT(py2.7)
species=Zebrafish
cpat.py -g nc.fa -d ~/biosoft/CPAT-1.2.4/dat/${species}_logitModel.RData -x ~/biosoft/CPAT-1.2.4/dat/${species}_Hexamer.tsv -o CPAT.out
cat CPAT.out|awk 'BEGIN{OFS="\t"}NR>1&&$6<=0.38{print $1,$6}'> CPAT_nc.list 

#合并
mv *nc.list result/
#提取gene list
#不同文库的情况
#cat CPAT_nc.list|awk '{split($1,x,"_");print tolower(x[1]) "_" x[2] "_" tolower(x[3])}' >CPAT.list

cat CPAT_nc.list|awk '{print tolower($1)}' > CPAT.list
cat PLEK_nc.list |awk '{split($3,x," ");print x[1]}'|sed 's/>//g' >PLEK.list
cat CNCI_nc.list |awk 'BEGIN{FS=" "}{print $1}'> CNCI.list

#取交集
comm -12 <( sort CNCI.list ) <( sort PLEK.list ) >CNCI_PLEK.list 
comm -12 <( sort CNCI_PLEK.list ) <( sort CPAT.list ) >CNCI_PLEK_CPAT.list 




##去掉比对上pfam库的序列
ln -s ../blast/result/pfam.anno.out
cat pfam.anno.out|cut -f 1 > pfam.list
#仅在CNCI_PLEK_CPAT.list 中出现的序列
comm -23 <( sort CNCI_PLEK_CPAT.list ) <( sort pfam.list ) >CNCI_PLEK_CPAT_pfam.list 

#提取lncRNA序列
cat ../tofu.collapsed.rep.fa |seqkit grep -f CNCI_PLEK_CPAT.list > lncRNA.original.fasta

#ORF预测 有编码潜能的序列
dumb_predict.py lncRNA.original.fasta angel.dumb --cpus 15
#--min_aa_length 300 
cat angel.dumb.final.cds|grep ">" |cut -d " " -f 1|sed 's/>//g'|sed 's/|.*//g' > angel.list
#去掉angel.list中的序列(23条)
cat lncRNA.original.fasta |seqkit grep -f angel/angel.list -v > lncRNA.fasta
##
#ORFfinder
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
~/biosoft/ORFfinder -in polish_assembly.fasta


#新转录本 262-23=239
cat lncRNA.fasta|grep ">" |cut -d " " -f 1 |sed 's/>//g' > lncRNA.final.list
#取差集 840-239=601
comm -23 <( sort ../no_genes.list ) <( sort lncRNA.final.list ) >novel.tr.list 
#提取序列
cat ../tofu.collapsed.rep.fa |seqkit grep -f novel.tr.list > novel.tr.fasta






#############transdecoder
mkdir Transdecoder && cd Transdecoder
ln -s ../lncRNA/result/novel.tr.fasta .
TransDecoder.LongOrfs -t novel.tr.fasta
TransDecoder.Predict -t novel.tr.fasta

cat novel.tr.fasta.transdecoder.cds|awk '
BEGIN{FS=OFS=" "}
/>/{split($7,x,":");
s1=">" x[1]":"x[2];
s2=$4;
s3=$5;
gsub("[\\(\\)]","",$6);
s4="strand:" $6;
split(x[3],y,"(");
s5="pos:" y[1]
print s1,s2,s3,s4,s5
}
!/>/{print}' > novel.cds.fa

cat novel.tr.fasta.transdecoder.pep|awk '
BEGIN{FS=OFS=" "}
/>/{split($7,x,":");
s1=">" x[1]":"x[2];
s2=$4;
s3=$5;
gsub("[\\(\\)]","",$6);
s4="strand:" $6;
split(x[3],y,"(");
s5="pos:" y[1]
print s1,s2,s3,s4,s5
}
!/>/{print}' > novel.pep.fa




mkdir bestHits
mkdir out
mkdir result
#############blast
#转录本功能注释是使用blastx，转录本是cDNA，跟蛋白数据库比对
#nr,trembl,pfam
#sp
ln -s ../lncRNA/result/novel.tr.fasta .

cd analysis/pro14_Maize/blast/
cDNA=novel.tr.fasta
db=sp #change

dbdir=${db}dir
dbdir=${!dbdir}
log=$db.log
out=$db.out
echo $dbdir $cDNA $log $out

nohup blastx -db $dbdir -query $cDNA -num_threads 20 -out $out -outfmt "6" -evalue 1e-5 -max_target_seqs 1 > $log 2>&1 &






#############SSR
#MISA 输出到fasta所在目录：nc.fa 
MISADIR=/home/liyb/biosoft/MISA
cp $MISADIR/misa.ini ./
ln -s ../lncRNA/nc.fa
perl $MISADIR/misa.pl nc.fa
mv nc.fa.misa SSR.misa
mv nc.fa.statistics SSR.statistics
