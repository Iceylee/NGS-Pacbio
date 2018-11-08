mkdir repeatmasker
mkdir trf
mkdir sRNA
mkdir tRNA
mkdir rRNA
mkdir islands
mkdir prophage
mkdir transposon

#############1.Repeat（repeatmaskter & trf）#############
# 物种找不到的话，不写物种参数。
cd repeatmasker
RepeatMasker -parallel 12 -html -gff -dir ./ ../assembly.fasta 
# 添加repeat class
sed '1,3 d' assembly.fasta.out|awk 'BEGIN{FS=" ";OFS="\t"}{locus=$5"-"$6"-"$7;print locus,$11}' > locus_repeatfamily.tab
awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}!/#/{locus=$1"-"$4"-"$5;print $0,a[locus]}' locus_repeatfamily.tab assembly.fasta.out.gff |grep -v "#" |awk 'NR==1{print "Sequence\tSource\tMethod\tStart\tEnd\tScore\tStrand\tPhase\tTarget\tRepeat_class"}NR>1{print $0}' >RepeatMasker.txt


# 去掉 |quiver
# cat temp.gff|awk 'BEGIN{OFS="\t"}{split($1,x,"|");print x[1],$2,$3,$4,$5,$6,$7,$8,$9,$10}'> temp2.gff


#最终结果不只包含散在重复序列。其实大量都是简单重复。
cat repeatmasker.out |cut -f 10|sort -k1,1 -u #种类

bedtools getfasta -fi ../assembly.fasta -bed assembly.fasta.out.gff > repeatmasker.fasta



# trf
cd trf
trf ../assembly.fasta 2 7 7 80 10 50 500 -f -d -m
mv *.dat out.dat
#将文件按分隔符分成多个
awk '{print $0 "Sequence: "> "file" NR}' RS='Sequence: ' out.dat
#选择大于4k的file：summary文件确认有repeat的contig数目 4
#将包含unitig的行和起始为数字的行 ，传给awk。第一行的unitig信息存为第一列，后面各行添加第一列信息。
for file in `ls file*`
do
grep -E '^[0-9]|unitig' $file|awk 'BEGIN{FS=" ";OFS="\t";}NR==1{name=$1}{$1=$1}NR>1{print name"\t"$0}' > $file.out
done
cat *.out >> all.out
#check with summary
#wc -l trf.out #1893
#去掉ACGT四列
cat all.out|cut -f 1,2,3,4,5,6,7,8,9,14,15,16|sed "1 i Sequence\tStart\tEnd\tPeriodSize\tCopyNumber\tConsensusSize\tPercentMatches\tPercentIndels\tScore\tEntropy\tConsensusSequences\tRepeatSequences" > Trf.txt








#############2.ncRNA#############
###3.1 sRNA(Rfam+cmscan)
	#参数Z = 11M *2
#耗时需要等
cd sRNA
nohup cmscan -Z 22 --cut_ga --rfam --nohmmonly --tblout my-genome.tblout --fmt 2 --cpu 20 --clanin ~/biosoft/Rfam/Rfam12.2.claninfo ~/biosoft/Rfam/Rfam.cm ../assembly.fasta > my-genome.cmscan &


grep -v "#" my-genome.tblout|awk '
BEGIN {FS=" "; OFS="\t"}
{if($20 != "=")print $2,$3,$4,$10,$11,$12,$17,$18}'|sed '1 i target_name\taccession\tquery_name\tquery_start\tquery_end\tstrand\tscore\tEvalue' > my-genome.tblout.final.xls

#统计各RNA类别个数
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type=="") type="Others"; count[type]+=1;}END{for(type in count) print type, count[type];}' ~/analysis/database/rfam_anno_class.txt my-genome.tblout.final.xls > ncRNA.stats

#统计各RNA类别总长度
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type=="") type="Others"; if($6=="-")len[type]+=$4-$5+1;if($6=="+")len[type]+=$5-$4+1}END{for(type in len) print type, len[type];}' ~/analysis/database/rfam_anno_class.txt my-genome.tblout.final.xls > ncRNA.length

#输出sRNA的全部记录
awk 'BEGIN{OFS=FS="\t"}ARGIND==1{a[$1]=$4;}ARGIND==2{type=a[$2]; if(type==" sRNA") print $0;}' ~/analysis/database/rfam_anno_class.txt my-genome.tblout.final.xls > sRNA.out
#加入第二列feature：sRNA
cat sRNA.out| awk 'BEGIN{FS=OFS="\t";}{print $3,"sRNA",$4,$5,$6,$1; }'|sed '1 i #Sequence\tFeature\tStart\tEnd\tStrand\tTarget_name'  >sRNA.gff
#调整至gff格式。并将start end顺序调整，end永远>start
cat sRNA.gff|awk '
BEGIN{FS=OFS="\t"}
FNR!=1{
	if ($5=="+") {start=$3;end=$4}
	if ($5=="-") {start=$4;end=$3}
	print $1,"rfam",$2,start,end,"-",$5,$6}
FNR==1{print $1,"Source",$2,$3,$4,"Score",$5,$6}

' > sRNA.gff 

bedtools getfasta -fi ../assembly.fasta -bed sRNA.gff > sRNA.fasta


#linxia
perl /bin/rfam_scan.pl -blastdb /biosoft/rfam11.0/Rfam.fasta /biosoft/rfam11.0/Rfam.cm polish_assembly.fasta -o polish_assembly.rfam.gff3 
python sumNCRNA.py -r polish_assembly.rRNA.gff2 -t polish_assembly.tRNA.out -s polish_assembly.rfam.gff3 -o polish_assembly.ncRNA.result.csv 

###3.2 tRNA(tRNAscan)
cd tRNA
tRNAscan-SE --bact -o tRNA.out -f tRNA.ss -m tRNA.stats ../assembly.fasta
cat tRNA.ss |grep Length | awk '{split($2,x," ");num = x[2];len += num}END{print len}'
     #计算总长度 !4566
cat tRNA.ss|awk '/unitig_0/{split($0,x,":");split(x[2],y," ");sum=sum+y[1]}END{print sum}'|head
cat tRNA.out|awk '
BEGIN{FS=OFS="\t"}
NR>3{print $1,"tRNAscan-SE","tRNA",$3,$4,$9,"+",$5,$6}'|sed "1 i #Sequence\tSource\tFeature\tStart\tEnd\tScore\tStrand\tType\tCodon" > tRNA.gff
#不知为何报错
bedtools getfasta -fi ../assembly.fasta -bed tRNA.gff > tRNA.fasta


cat tRNA.ss|awk '
/Seq/{
gsub(/Seq: /,"",$0)
print $0
}
/unitig/{
print ">"$0
}' > tRNA.fasta


###3.3 rRNA(rRNAmmer)
#-S arc bac euk
#10min
cd rRNA
#perl /home/liyb/biosoft/rRNAmmer/rnammer -S bac -m lsu,ssu,tsu -gff - < ../assembly.fasta > out.gff
perl /home/liyb/biosoft/rRNAmmer/rnammer -S bac -multi -f rRNA.fna -h rRNA.hmmreport -xml rRNA.xml -gff out.gff ../../assembly.fasta

grep -v "#" out.gff|awk 'BEGIN{OFS="\t"}{print $1,"rRNAmmer",$3,$4,$5,$6,$7,$8,$9}' |sed "1 i #Sequence\tSource\tFeature\tStart\tEnd\tScore\tStrand\tFrame\tAttribute" > rRNA.gff

bedtools getfasta -fi ../assembly.fasta -bed rRNA.gff > rRNA.fasta


#############3.islands（dimob）#############
cd islands

#gbk文件不能带多余的点号，后缀gbk
#多个contig会有报错信息，运行完为空文件，将gff和fasta按contig分开
for i in 0
do
cat ../assembly.fasta|seqkit grep -p unitig_${i} > unitig_${i}.fasta
cat ../prodigal/Thermoactinomyces_sp.CDF.gff|grep unitig_${i} > unitig_${i}.gff
python ~/analysis/scripts/gff_to_genbank.py unitig_${i}.gff unitig_${i}.fasta
mv unitig_${i}.gb unitig_${i}.gbk
~/biosoft/islandpath-master/Dimob.pl unitig_${i}.gbk ${i}_GIs.txt
done

for file in `ls *_GIs.txt`
do
	name=${file/_GIs.txt/}
cat ${file} | awk -v contig="$name" 'BEGIN{OFS="\t"}{print "unitig_"contig,$0}' >> all.GIs.out
done


cat all.GIs.out|awk 'BEGIN{OFS="\t"}{$2="Genomic Island";print $0}' |sed "1 i Sequence\tFeature\tStart\tEnd" > islands.txt

# gff提取gene坐标
ln -s ../prodigal/Thermoactinomyces_sp.CDF.gff gff
awk 'BEGIN{OFS="\t"}
NR==FNR&&FNR>1{i++;sequence[i]=$1;start[i]=$3;end[i]=$4;next}
{for (i in start) 
	{if($4>=start[i]&&$5<=end[i]) print sequence[i],"Dimob","Gene",$4,$5,".",$7,$9,start[i],end[i],end[i]-start[i]+1}}' islands.txt gff> temp
cat temp | awk 'BEGIN {OFS="\t"}
	{print $0}'|sed "1 i #Sequence\tSource\tFeature\tGene_Start\tGene_End\tScore\tStrand\tGene_ID\tIsland_Start\tIsland_End\tIsland_Length" > GeneIslands.txt

bedtools getfasta -fi ../assembly.fasta -bed GeneIslands.txt > GeneIslands.fasta

#############4.transposon（TransposonPSI）#############
#若蛋白序列 最后一个参数 prot
cd transposon
~/biosoft/TransposonPSI_08222010/transposonPSI.pl ../assembly.fasta nuc
cat assembly.fasta.TPSI.allHits.chains.gff3|awk '{len=$5-$4;sum=sum+len}END{print sum}' #!39622

sed "1 i Sequence\tSource\tMethod\tStart\tEnd\tScore\tStrand\tPhase\tAttributes" assembly.fasta.TPSI.allHits.chains.gff3 >Transposon.txt

bedtools getfasta -fi ../assembly.fasta -bed assembly.fasta.TPSI.allHits.chains.gff3 > Transposon.fasta


#######5.prophage#####
#在线phaser
cd prophage
#仅unitig0有结果
cat summary.txt|sed '1,31 d'|grep -v "\-----"|awk 'BEGIN{FS=" ";OFS="\t"}{sub(/[[:blank:]]+/,"",$0);print $0}' > unitig_0.txt
cat unitig_0.txt|sed "1d"|awk 'BEGIN{FS=" ";OFS="\t"}{$1=$1}1'|awk 'BEGIN{OFS="\t"}{print "unitig_0",$0}' > out
cat out |awk '
BEGIN{OFS="\t"}
	{$3="Prophage";
	$2="PHASTER";

	split($6,x,"-");
	print $1,$2,$3,x[1],x[2],".","+",$3,$4,$15}' >temp.out

cat temp.out|awk '
BEGIN{
	print "#Sequnece\tSource\tFeature\tStart\tEnd\tScore\tStrand\tRegion_Length\tCompleteness(score)\tMost_Common_Phage_Name(hit_genes_count)\t"}1
' > Prophage.txt

bedtools getfasta -fi ../assembly.fasta -bed Prophage.txt > Prophage.fasta


####crispr
cat CRISPR.txt|awk '
BEGIN{FS=OFS="\t"}
FNR!=1{

	print $1,"CRISPR","CRISPR_Finder",$3,$4,".","+",$5,$6}
FNR==1{print "#"$1,"Feature","Source",$3,$4,"Score","Strand",$5,$6}

' > CRISPR.gff 

bedtools getfasta -fi ../assembly.fasta -bed CRISPR.gff > CRISPR.fasta

####6.modification
gunzip *.gz
mv *modifications.gff modifications.gff
mv *modifications.csv modifications.csv
mv *motifs.gff motifs.gff
mv *motif_summary.csv motif_summary.csv

