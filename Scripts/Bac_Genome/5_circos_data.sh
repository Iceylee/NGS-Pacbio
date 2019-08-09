# assembly.fasta
# rRNA.gff
# tRNA.gff
# sRNA.gff
# *.gff
# Prophage.txt
# CRISPR.txt
# islands.txt

#0.karyotype
#得到长度，去掉第一行；按数字排序；匹配chr来设置颜色
seqkit fx2tab -l -g -n -i -H assembly.fasta |cut -f 1,4|sed '1 d'|sort -V -k1,1|awk '{i=i+1;color="chr"i;print "chr - "$1" "$1" 0 "$2" "color}' > karyotype.txt

GC=??
#1.GC content & GC skew
python ~/work/scripts/GCcalc.py -f assembly.fasta -w 10000 -s 10000 >gc_content_skew.txt
##############change 参数
cat gc_content_skew.txt|awk -v var="$GC" '{OFS="\t"}{print $1,$2,$3,$4-var}'>gc_count.txt
cat gc_content_skew.txt|awk '{OFS="\t"}{print $1,$2,$3,$5}'>gc_skew.txt


#2.CDS:正负义链
gff=Thermoactinomyces_sp.CDF.gff
cat $gff |awk '$7=="+"&&FNR>1{OFS="\t";print $1,$4,$5}' >sensegene.gff
cat $gff |awk '$7=="-"&&FNR>1{OFS="\t";print $1,$4,$5}' >antigene.gff

#color：color.list即希望着色的颜色列表 
#_a4 为透明度
awk 'NR==FNR{i=i+1;a[i]=$1;next}{num=int(rand()*7)+1;print $1,$2,$3,"fill_color="a[num]"_a2"}' ~/work/scripts/color2.list sensegene.gff >sensegene_color.gff
awk 'NR==FNR{i=i+1;a[i]=$1;next}{num=int(rand()*7)+8;print $1,$2,$3,"fill_color="a[num]"_a2"}' ~/work/scripts/color2.list antigene.gff >antigene_color.gff

#3.ncRNA.gff
# cat *RNA.gff|grep -v "Sequence"|cut -f 1,2,3,4 >nc.gff
cat rRNA.gff|grep -v "Sequence"|cut -f 1,4,5 >rRNA_plot.gff
cat tRNA.gff|grep -v "Sequence"|cut -f 1,4,5 >tRNA_plot.gff
cat sRNA.gff|grep -v "Sequence"|cut -f 1,4,5 >sRNA_plot.gff

#4.
cat Prophage.txt|sed '1d'|cut -f 1,3,4 > Prophage_plot.txt
cat CRISPR.txt|sed '1d'|cut -f 1,3,4 > CRISPR_plot.txt
cat islands.txt|sed '1d'|cut -f 1,3,4 > GIs_plot.txt

#5.文件夹
circos -conf plot/circos.conf