genome=$1 #from arguments

echo "##########stat##########" > report.txt
seqkit stat assembly.fasta >report.txt
echo "##########gc##########" >> report.txt
cat prodigal/prodigal.stat|grep gc_cont >> report.txt

#基因个数
num="$(cat prodigal/prodigal.cds|grep ">"|wc -l )"
# gene总长  
len="$(grep -v "#" prodigal/prodigal.gff|awk '$5>$4{len=len+$5-$4+1}$5<$4{len=len+$4-$5+1}END{print len}')"
echo "##########prodigal##########" >> report.txt
echo -n "gene number: " >> report.txt
echo $num >>report.txt
echo -n "gene length sum: " >> report.txt
echo $len >>report.txt
echo -n "gene average len: " >>report.txt
echo $len/$num|bc -l >>report.txt
echo -n "gene percentage: " >>report.txt
echo $len/$genome|bc -l >> report.txt


# repeatmasker:out文件不只散在重复。
echo "##########repeatmasker##########" >> report.txt
cat repeatmasker/assembly.fasta.tbl|grep "SINEs" >> report.txt
cat repeatmasker/assembly.fasta.tbl|grep "LINEs" >> report.txt
cat repeatmasker/assembly.fasta.tbl|grep "LTR elements" >> report.txt
cat repeatmasker/assembly.fasta.tbl|grep "DNA elements" >> report.txt
cat repeatmasker/assembly.fasta.tbl|grep "Unclassified" >> report.txt

cat repeatmasker/assembly.fasta.tbl|grep "Total interspersed repeats" >> report.txt


# trf
# number
num="$(sed '1 d' trf/Trf.txt|wc -l)"
# total Length
len="$(sed '1 d' trf/Trf.txt|awk '{len=$3-$2+1+len}END{print len}')"
echo "##########trf##########" >> report.txt
echo -n "trf number: " >> report.txt
echo $num >>report.txt
echo -n "trf length sum: " >> report.txt
echo $len >>report.txt
echo -n "genome percentage: " >>report.txt
echo $len/$genome|bc -l >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt

# rRNA
# number
num="$(sed '1 d' rRNA/rRNA.gff |wc -l)"
# total Length
len="$(sed '1 d' rRNA/rRNA.gff|awk '{len=len+$4-$3+1}END{print len}')"
echo "##########rRNA##########" >> report.txt
echo -n "rRNA number: " >> report.txt
echo $num >>report.txt
echo -n "rRNA length sum: " >> report.txt
echo $len >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt


#sRNA
echo "##########sRNA##########" >> report.txt
echo -n "sRNA number: " >> report.txt
cat sRNA/ncRNA.stats|grep sRNA >>report.txt
echo -n "sRNA length sum: " >> report.txt
cat sRNA/ncRNA.length|grep sRNA>>report.txt


#tRNA
# number
num="$(sed "1 d" tRNA/tRNA.gff|wc -l)"
# total Length
len="$(cat tRNA/tRNA.ss |grep Length | awk '{split($2,x," ");num = x[2];len += num}END{print len}')"
echo "##########tRNA##########" >> report.txt
echo -n "tRNA number: " >> report.txt
echo $num >>report.txt
echo -n "tRNA length sum: " >> report.txt
echo $len >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt

#islands
# number
num="$(cat islands/islands.txt |sed '1 d' |wc -l )"
# total Length
len="$(cat islands/islands.txt|awk '$4>$3{sum=sum+$4-$3}$4<$3{sum=sum+$3-$4}END{print sum}')"
echo "##########islands##########" >> report.txt
echo -n "islands number: " >> report.txt
echo $num >>report.txt
echo -n "islands length sum: " >> report.txt
echo $len >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt

#transposon
# number
num="$(cat transposon/Transposon.txt |sed '1 d' |wc -l )"
# total Length
len="$(sed '1 d' transposon/Transposon.txt|awk '$5>$4{len=len+$5-$4+1}$5<$4{len=len-$5+$4+1}END{print len}')"
echo "##########transposon##########" >> report.txt
echo -n "transposon number: " >> report.txt
echo $num >>report.txt
echo -n "transposon length sum: " >> report.txt
echo $len >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt


#prophage
# number
num="$(sed '1 d' prophage/Prophage.txt|wc -l)"
# total Length
len="$(sed '1 d' prophage/Prophage.txt |awk '$4>$3{len=len+$4-$3+1}$4<$3{len=len-$4+$3+1}END{print len}')"
echo "##########prophage##########" >> report.txt
echo -n "number: " >> report.txt
echo $num >>report.txt
echo -n "length sum: " >> report.txt
echo $len >>report.txt
echo -n "average len: " >>report.txt
echo $len/$num|bc -l >> report.txt


#CRISPR
# number
# num="$(sed '1 d' CRISPR/CRISPR.txt|wc -l)"
# # total Length
# len="$(sed '1 d' CRISPR/CRISPR.txt|awk '$4>$3{len=len+$4-$3+1}$4<$3{len=len-$4+$3+1}END{print len}')"
# echo "##########CRISPR##########" >> report.txt
# echo -n "number: " >> report.txt
# echo $num >>report.txt
# echo -n "length sum: " >> report.txt
# echo $len >>report.txt
# echo -n "average len: " >>report.txt
# echo $len/$num|bc -l >> report.txt

######甲基化
#number
echo "##########modifications##########" >> report.txt
#求平均coverage
echo -n "coverage&num m4C: " >> report.txt
grep -v "#" modification/modifications.gff|grep "m4C"|awk '{match($0,/coverage=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt
echo -n "coverage&num m6A: " >> report.txt
grep -v "#" modification/modifications.gff|grep "m6A"|awk '{match($0,/coverage=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt
echo -n "coverage&num modified_base: " >> report.txt
grep -v "#" modification/modifications.gff|grep "modified_base"|awk '{match($0,/coverage=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt


#求平均QV
echo -n "QV m4C: " >> report.txt
grep -v "#" modification/modifications.gff|grep "m4C"|awk '{match($0,/identificationQv=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt
echo -n "QV m6A: " >> report.txt
grep -v "#" modification/modifications.gff|grep "m6A"|awk '{match($0,/identificationQv=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt
echo -n "QV modified_base: " >> report.txt
grep -v "#" modification/modifications.gff|grep "modified_base"|awk '{match($0,/identificationQv=[0-9]+/,m);match(m[0],/[0-9]+/,n);sum=sum+n[0];i=i+1}END{print i,sum/i}' >> report.txt