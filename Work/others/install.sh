#到https://sourceforge.net/projects/tgicl/files/?source=navbar下载最新版的TGICL (以最新版v2.1版为例)。建议下载TGICL-2.1.tar.gz原程序

tar -zxvf TGICL-2.1.tar.gz

cd TGICL-2.1

perl Build.PL
./Build
./Build test
./Build install

#如果哪一步有permission denied的提示，则需要在命令前加sudo

#帮助信息
perldoc -F /data1/software/miniconda2/bin/tgicl

tgicl -F ../Trinity.fasta -c 15
#Use of :locked is deprecated at /data1/software/miniconda2/lib/perl5/site_perl/5.22.0/TGI/DBDrv.pm line 36. 该错误可忽略

#error
sh: /data1/software/miniconda2/bin/psx: /lib/ld-linux.so.2: bad ELF interpreter: No such file or directory
#据说是64位系统安装32位软件导致
sudo yum install glibc.i686
sudo yum install zlib.i686

