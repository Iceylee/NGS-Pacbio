http://bioinf.uni-greifswald.de/augustus/binaries/


# maker
less INSTALL
cd maker/src/
perl Build.PL 
# 查看需要安装的perl模块 部分软件还未安装，因此询问MPI时候选择N
# src下locations为各软件的安装路径。也可手动下载
./Build status #查看目前缺少哪些
#安装 perl 模块
./Build installdeps #询问perl 多线程是否指向foks，选择yes
#安装依赖程序
./Build exonerate
./Build mpich2 #mpi软件，最好是管理员安装。
#mpich2后
perl Build.PL  #yes。默认路径



#安装maker
./Build install

#安装完成 出现bin目录
#添加入环境变量






#genemaker-es
gunzip gm_key_64.gz
mv gm_key_64 ~/.gm_key
less README.GeneMark-ES-suite

#Required PERL modules available on CPAN
   # YAML
   # Hash::Merge
   # Logger::Simple
   # Parallel::ForkManager
#安装perl模块



#RAPSearch
./install


#GeneMark-ES/ET
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Qg87n/gm_et_linux_64.tar.gz
tar xf gm_et_linux_64.tar.gz
mv gm_et_linux_64/gmes_petap/ ~/opt/biosoft
wget http://topaz.gatech.edu/GeneMark/tmp/GMtool_Qg87n/gm_key_64.gz
gzip -dc gm_key_64.gz > ~/.gm_key
cpan YAML Hash::Merge Logger::Simple Parallel::ForkManager
#模块有点问题 安装了但还是无法找到 最后加入了软件的lib路径/data1/software/quast-4.6.0/quast_libs/genemark-es/lib/ 在 $PERL5LIB

#EVM
cd ~/src
wget -4 https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz
tar xf v1.1.1.tar.gz
mv EVidenceModeler-1.1.1 ~/opt/biosoft/

#augustus
source activate annotation
conda install augustus=3.3


#exonerate
cd ~/src
wget http://ftp.ebi.ac.uk/pub/software/vertebrategenomics/exonerate/exonerate-2.2.0-x86_64.tar.gz
tar xf exonerate-2.2.0-x86_64.tar.gz
mv exonerate-2.2.0-x86_64 ~/opt/biosoft/exonerate-2.2.0# .bashrc添加环境变量export PATH=~/opt/biosoft/exonerate-2.2.0:$PATH# 或conda install -c bioconda exonerate


#maker
conda create -n marker marker

#BRAKER
cpan File::Spec::Functions Module::Load::Conditional POSIX Scalar::Util::Numeric YAML File::Which Logger::Simple Hash::Merge
wget -4 http://exon.biology.gatech.edu/GeneMark/Braker/BRAKER2.tar.gz
tar xf BRAKER2.tar.gz -C ~/opt/biosoft

echo "export PATH=$PATH:$HOME/opt/biosoft/BRAKER_v2.1.0/" >> ~/.bashrc
# 在~/.bashrc设置如下软件所在环境变量
export AUGUSTUS_CONFIG_PATH=/data1/software/miniconda2/envs/annotation/config/
export AUGUSTUS_BIN_PATH=/data1/software/miniconda2/envs/annotation/bin
export AUGUSTUS_SCRIPTS_PATH=/data1/software/miniconda2/envs/annotation/bin/
export BAMTOOLS_PATH=/data1/software/miniconda2/envs/annotation/bin/
export GENEMARK_PATH=/data1/software/gm_et_linux_64/gmes_petap/
export SAMTOOLS_PATH=/usr/local/bin/
export BLAST_PATH=/data1/software/miniconda2/bin/



##RepeatModel
# RECON
cd ~/src
wget -4 http://repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
tar xf RECON-1.08.tar.gzcd RECON-1.08/src
make && make installcd ~/src
mv RECON-1.08 ~/opt/biosoft
# nesg
cd ~/src
mkdir nesg && cd nesg
wget -4 ftp://ftp.ncbi.nih.gov/pub/seg/nseg/*
make
mv nmerge nseg ~/opt/bin/
# RepeatScout
http://www.repeatmasker.org/RepeatScout-1.0.5.tar.gz
# RepeatModel
wget -4 http://repeatmasker.org/RepeatModeler/RepeatModeler-open-1.0.11.tar.gz
tar xf RepeatModeler-open-1.0.11.tar.gz
mv RepeatModeler-open-1.0.11 ~/opt/biosoft/cd ~/opt/biosoft/RepeatModeler-open-1.0.11
# 配置
perl ./configure
export PATH=~/opt/biosoft/maker:$PATH


##GenomeThreader
wget -4 http://genomethreader.org/distributions/gth-1.7.0-Linux_x86_64-64bit.tar.gz
tar zxvf gth-1.7.0-Linux_x86_64-64bit.tar.gz
# 修改.bashrc增加如下行
export PATH=$PATH:$HOME/opt/biosoft/gth-1.7.0-Linux_x86_64-64bit/bin
export BSSMDIR="/data1/software/gth-1.7.0-Linux_x86_64-64bit/bin/bssm"
export GTHATADIR="/data1/software/gth-1.7.0-Linux_x86_64-64bit/bin/gthdata"

#PASA
cpan DB_File URI::Escape DBI DBD::SQLite
# GMAP
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2017-11-15.tar.gz
tar xf gmap-gsnap-2017-11-15.tar.gzcd gmap-2017-11-15./configure --prefix=$HOME/opt/gmap
make && make install
# BLAT
cd ~/opt/bin
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && chmod 755 ./blat# Fasta3wget -4 http://faculty.virginia.edu/wrpearson/fasta/fasta36/fasta-36.3.8g.tar.gz && \
tar zxvf fasta-36.3.8g.tar.gz && \cd ./fasta-36.3.8g/src && \
make -f ../make/Makefile.linux_sse2 all && \
cp ../bin/fasta36 ~/opt/bin
# 以上程序需添加到环境变量中
# PASApipeline
git clone https://github.com/PASApipeline/PASApipeline.git
cd PASApipeline 
git checkout feature/sqlite 
git submodule init 
git submodule update 
make

#gmap
wget http://research-pub.gene.com/gmap/src/gmap-gsnap-2018-05-30.tar.gz
tar zvxf gmap-gsnap-2018-05-30.tar.gz
cd gmap-2018-05-30
./configure 
##--prefix=$HOME/opt/gmap
sudo make && sudo make install

##genewise
wget http://www.ebi.ac.uk/~birney/wise2/wise2.4.1.tar.gz
tar zxfv wise2.4.1.tar.gz


#yum install *glib*
sudo yum install glib
sudo yum install glib-devel

perl -p -i -e 's/getline/get_line/g' ./HMMer2/sqio.c
perl -p -i -e 's/isnumber/isdigit/' models/phasemodel.c

#报错：glib-config command not found则 把src下面所有makefile中的glib-config字段都替换为pkg-config --libs glib-2.0
#find ./ -name makefile | xargs sed -i 's/glib-config/pkg-config --libs glib-2.0/'

#export C_INCLUDE_PATH=/usr/include/glib-2.0/:/usr/lib64/glib-2.0/include/:$C_INCLUDE_PATH


make all
export WISECONFIGDIR=/data1/software/wise2.4.1/wisecfg/
make test
#有报错 但test passed就行


echo 'PATH=$PATH:/data1/software//wise2.4.1/src/bin/' >> ~/.bashrc
echo 'export WISECONFIGDIR=/data1/software/wise2.4.1/wisecfg/' >> ~/.bashrc 
source ~/.bashrc






