#orthomcl pipeline
git clone https://github.com/apetkau/orthomcl-pipeline.git

#手动安装
perl -MCPAN -e shell

#install cpanm
#https://perlbrew.pl/Perlbrew-中文簡介.html
perlbrew install-cpanm
cpanm BioPerl DBD::mysql DBI Parallel::ForkManager YAML::Tiny Set::Scalar Text::Table Exception::Class Test::Most Test::Warn Test::Exception Test::Deep Moose SVG Algorithm::Combinatorics

conda install perl-bioperl perl-test-warn 



#分析
#sql 创建叫orthomcl的database
create user 'orthomcl'@localhost identified by 'BYzh-123';
create database orthomcl;
grant all on orthomcl.* to 'orthomcl'@'localhost';
flush privileges;
show databases;


#准备mysql配置文件
#cp /data1/software/miniconda2/share/orthomcl/orthomcl.config.template /home/liyubing/analysis/0528_Orthomcl/orthomcl.config
orthomcl-setup-database.pl --user orthomcl --password BYzh-123 --host localhost --database orthomcl --outfile configure_outfile.conf --no-create-database
#报错找不到/tmp/mysql.sock 解决办法
ln -s /private/var/mysql/mysql.sock /tmp/mysql.sock

for file in `ls *.fasta`
do
	echo "sed -i \"s/-/_/g\" $file"
done


#开始pipeline分析
orthomcl-pipeline -i annotations-small/ -o orthomcl-output-small -m configure_outfile.conf --nocompliant


#venn diagram
##准备genome-groups.txt
nepal: VC_25,VC_26,VC_14
haiti: 2010EL_1749,2010EL_1786,2010EL_1796,2010EL_1798,2011EL_2317,2012V_1001,3554_08,VC_10,VC_15,VC_18,VC_19,VC_1,VC_6
c6706: C6706

nml_parse_orthomcl.pl -i orthomcl-output-small/groups/groups.txt -g genome-groups.txt -s --draw -o orthomcl-stats.txt --genes





#报错原因
文件名不能带 破折号。
后续取名字的时候，没有取完整。

#替换文件名 破折号替换成下划线
for f in *.fasta
do 
	mv -- "$f" "${f//-/_}"; 
done











#islands
cpanm Data::Dumper Log:Log4perl Config::Simple Moose MooseX::Singleton 
#Bio::Perl 手动下载
cd BioPerl-1.6.0
perl Build.PL  --install_base /sam/bioperl
./Build test
./Build install
#gff_to_gbk.py
pip install bcbio-gff


##kobas
conda info --envs
source activate kobasPy2.7
conda install sqlite biopython numpy pandas rpy2 matplotlib

R
source("https://bioconductor.org/biocLite.R")
biocLite("qvalue")
biocLite("GSVA")
biocLite("gage")
biocLite("EnrichmentBrowser")
q()

cd ~/biosoft/kobas/
R CMD INSTALL EnrichmentBrowser_2.4.5_CA_edit.tar.gz

R
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")
biocLite(c("ComplexHeatmap","ReportingTools","Rgraphviz","biocGraph"))


install.packages("kernlab")
source("https://bioconductor.org/biocLite.R")
biocLite("RBGL")
biocLite("Rgraphviz")
install.packages("flexmix")
install.packages("fpc")
conda install -c r r-xml=3.98_1.5
biocLite("Category")
biocLite("pathview")
biocLite("biocGraph")
install.packages("dendextend")
biocLite("GOstats")

biocLite("ggbio")
biocLite("ComplexHeatmap")
biocLite("ReportingTools")
biocLite("EnrichmentBrowser")
install.packages("RSQLite")
biocLite("OrganismDbi")


wget https://pypi.python.org/packages/42/02/981b6703e3c83c5b25a829c6e77aad059f9481b0bbacb47e6e8ca12bd731/pysqlite-2.8.3.tar.gz#md5=033f17b8644577715aee55e8832ac9fc
tar -zxf pysqlite-2.8.3.tar.gz
cd pysqlite-2.8.3
#edit the setup.cfg file: change include_dirs and library_dirs to the PATH you want to install.)
python setup.py build
python setup.py install

##install islander(下载)
#将prodigal和hmmsearch的路径加入bin
#安装pftools(下载bin，直接可运行)
ftp://ftp.lausanne.isb-sib.ch/pub/software/unix/pftools/pft2.3/executables.
perl Islander.pl example/NC_000913 --verbose --reisland --table 11 --nocheck


