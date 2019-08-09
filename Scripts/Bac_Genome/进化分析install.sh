#CEGMA 
git clone https://github.com/KorfLab/CEGMA_v2.git
make

export CEGMA="/home/liyubing/biosoft/CEGMA_v2/"
#export CEGMATMP="path"  #default /tmp
export PERL5LIB="$PERL5LIB:$CEGMA/lib"

bin/cegma 

#geneid 安装
wget ftp://genome.crg.es/pub/software/geneid/geneid_v1.4.4.Jan_13_2011.tar.gz
make
src/geneid

#pseudo-finder
git clone https://github.com/filip-husnik/pseudo-finder.git

#install Python3 miniconda https://conda.io/miniconda.html
conda install -c conda-forge biopython
#conda install -c bioconda blast

conda install numpy reportlab plotly pandas -y

conda create --name py3 python=3

#需要提前调整转录起始位点到0（如果是一个闭合环）




#ORI finder
https://github.com/ShanSabri/oriFinder.git
python oriFinder.py ecoli_genome.fasta




#shuailab-pseudogenepipe
git clone https://github.com/ShiuLab/PseudogenePipeline.git
#需要fasta 36版本
wget http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz
cd src
make -f ../make/Makefile.linux_sse2 all
#path：/home/liyubing/biosoft/fasta-36.3.8g/bin


#install PAML
#http://abacus.gene.ucl.ac.uk/software/paml.html
tar xf paml4.9h.tgz
rm bin/*.exe 
cd src 
make -f Makefile 
ls -lF 
rm *.o 
mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin 
cd .. 
ls -lF bin 
bin/baseml 
bin/codeml 
bin/evolver 

#KaKs_Calculator:http://bigd.big.ac.cn/tools/kaks/download
cd src
#以下三文件 首行加include行
KaKs.cpp #include "string.h"
AXTConvertor.cpp #include "stdlib.h"
GY94.cpp #include "string.h"

make 

#clustalw
wget http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz

tar -zxvf clustalw-2.1.tar.gz
./configure
make
make install

src/clustalw2


#pal2nal.pl 支持PAML的小工具
git clone https://github.com/mrmckain/FASTKs.git


#codon2kaks.pl

# 1gene,2seqs.
# pal2nal output => kaks_calculator input
<>;<>;
for(1..3){
$a=<>;$a=~s/\s+$//;
last unless(defined $a || $a eq '');
$b=<>;$b=~s/\s+$//;
<>;
($n1,$seq1)=(split/\s+/,$a);
($n2,$seq2)=(split/\s+/,$b);
$line1.=$seq1;$line2.=$seq2;
}
print "$n1,$n2\n$line1\n$line2\n";
print "\n";




