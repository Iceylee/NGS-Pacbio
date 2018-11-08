##############1.combine all hq&lq,and get fasta file##############
mv *hq* hq.fq
mv *lq* lq.fq
seqret -sequence hq.fq -outseq hq.fa
seqret -sequence lq.fq -outseq lq.fa
cat hq.fa >> lq.fa
mv lq.fa isoform.fa

ref_fa=GCF_000005005.2_B73_RefGen_v4_genomic.fna
ref_name=Maize_B73

##gmap索引
#修改 输出路径和物种名 9：29
# fasta-to-gmap-reference $ref_fa gmap_index $ref_name --organism organism
# pbservice import-dataset --host http://192.168.1.150 --port 8243 gmap_index/C_idella/gmapreferenceset.xml

##############2.gmap比对参考基因组########
#1.参考基因组索引  #1.3G 40min  #1h 900M
gmap_build -d $ref_name ref/$ref_fa
#geneapps 
#gmap_build -d gmap_db ref/alter.C_idella_female_scaffolds.fasta 

#2.比对  -D 索引所在目录 -d db前缀同文件夹名字 
#-t 8 cu04 2h #-t 20 
#-d gmap_db
mkdir gmap && cd gmap
nohup gmap -D /home/liyb/miniconda3/share/$ref_name -d $ref_name -f samse -n 0 -t 10 ../genome/isoform.fa > gmap.sam 2> gmap.sam.log &
#3.sort
#samtools view -b DM.sam -o DM.bam
#samtools sort -O SAM -o DM.sorted.sam DM.bam
sort -k 3,3 -k 4,4n gmap.sam > gmap.sorted.sam



##############3.ToFu########
source activate anapy2.7.9
export VENV_TOFU=~/biosoft/VENV_TOFU
source $VENV_TOFU/bin/activate

#1.融合基因
cd Fusion
fusion_finder.py --input ../genome/isoform.fa -s ../gmap/gmap.sorted.sam -o fusion 
#--cluster_report_csv cluster_report.csv cDNA-primer版本没有这个参数。csv在SMRTLink job 里面得到。

#2.非冗余转录本
cd collapse2
collapse_isoforms_by_sam.py --input ../genome/isoform.fa -s ../gmap/gmap.sorted.sam --dun-merge-5-shorter -o tofu
source deactivate 

#3.去冗余后重新比对和sort 
cd ../gmap2
nohup gmap -D /home/liyb/miniconda3/share/C_idella -d C_idella -f samse -n 0 -t 16 ../collapse/tofu.collapsed.rep.fa > gmap.collapsed.sam 2> gmap.collapsed.sam.log &
sort -k 3,3 -k 4,4n gmap.collapsed.sam > gmap.collasped.sorted.sam



#整理文件
mv gmap.collasped.sorted.sam collapsed.sam
gzip collapsed.sam -c > collapsed.sam.gz
gzip gmap.sorted.sam -c > gmap.sorted.sam.gz







##############4.MatchAnnot(py2.7)########
# 0.1 gff to gtf
#cufflinks:gffread
#####gffread CID_femalev11_genemodel_edit.gff -T -o my.gtf
#??有少量序列报错

# 0.2 为避免matchanno报错，除去gtf中没有的染色体序列
ln -s ../ref/CID_femalev11_genemodel_edit.gtf my.gtf
ln -s ../gmap2/gmap.collasped.sorted.sam ./
cat my.gtf |grep -v "#"|cut -f 1|sort -k1,1 -u > chr_in_gtf.list #32811
#sam文件的@行都去掉了，应该不影响？
awk 'NR==FNR{a[$1]=$1;next}($3 in a){print $0}' chr_in_gtf.list gmap.collasped.sorted.sam > chr_filter.sam 

#gtf文件
#报错：tranName = re.search(regexTran, attrs).group(1)
#AttributeError: 'NoneType' object has no attribute 'group'
#表示需要transcriptName

cat my.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[1],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print $0 " transcript_name \""a[0] "\";"} ' >my2.gtf

#exon number



# 1.非冗余转录本比对结果与参考序列注释比较 
source activate Icey

#加参数--format alt 不报错“缺少exon信息”
matchAnnot.py --gtf my2.gtf --format alt chr_filter.sam > annotations.out


#gtf文件有多少基因
cat my2.gtf|awk 'BEGIN{FS=OFS="\t"}
{
	split($9,x,";"); 
	match(x[2],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	print a[0]}'|sort -k1,1 -u|wc -l
#49316

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
cat temp5|cut -f 1,4,5 | sed "1 i PacBio_transcript_id\tAnnoted_gene_id\tAnnoted_transcript_id"> matchAnnot_result.txt




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
#处理原始gtf文件:第9行取第三个；匹配引号；去除引号；判断trID是否已在isoseq list中；无则输出

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {trList[$1]=1;next}
($3=="transcript" || $3=="exon") {
	split($9,x,";"); 
	match(x[3],/\".*\"/,a);
	gsub(/"/,"",a[0]);
	if (!(a[0] in trList)) print $0}
' isoseq.tr.list ../my2.gtf >ref.spec.gtf






#整理文件4：merge.gtf
#isoseq.info.gtf + ref.spec.gtf
cat ref.spec.gtf isoseq.info.gtf > merge.gtf






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









