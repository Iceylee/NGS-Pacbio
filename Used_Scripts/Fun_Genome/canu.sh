conda install canu
conda install busco

#MECAT
git clone https://github.com/xiaochuanle/MECAT.git
cd MECAT
make  #ERROR

#falcon
docker pull billandreo/falcon_ctr3
docker run -v /home/liyubing/analysis/test:/mydata -it billandreo/falcon_ctr3 #do-assemble.sh
docker run -it -rm billandreo/falcon_ctr3 /bin/bash

#smrtlink
##docker (ok)
docker pull jarvice/app-pacbio-smrtlink
sudo service docker start
docker run -it -rm jarvice/app-pacbio-smrtlink /bin/bash
#/opt/pacbio/smrtlink/smrtcmds/bin
pbsmrtpipe show-templates



#c 断点续传 q 安静模式
wget -c -q ftp://ftp.sra.ebi.ac.uk/vol1/ERA111/ERA1116568/bam/pb.bam

# build index for convert
~/opt/biosoft/smrtlink/smrtcmds/bin/pbindex pb.bam &
# convert bam to fasta
~/opt/biosoft/smrtlink/smrtcmds/bin/bam2fasta -o pb pb.bam &























#smrtlink install
#Paste the following into /etc/security/limits.conf just before the line # End of file:
*         hard    nofile      500000
*         soft    nofile      500000
root      hard    nofile      500000
root      soft    nofile      500000

#重开终端
ulimit -Sn #确认设置是否ok

#
sudo service docker start
sudo usermod -a -G docker liyubing

docker build --ulimit nofile=500000:500000 https://github.com/caseywdunn/docker-smrtlink.git#master -t smrtlink:5.0.1.9585
#-t 镜像的名字及标签，通常 name:tag 或者 name 格式
#报错，手动下载

git clone https://github.com/caseywdunn/docker-smrtlink.git

cd /home/liyubing/database/docker-smrtlink/docker
#在dockerfile所在目录执行
docker build --ulimit nofile=500000:500000 -t smrtlink:5.0.1.9585 .



#自行安装
# SMRT_USER=smrtanalysis
# sudo useradd -ms /bin/bash $SMRT_USER #添加用户
# usermod  -g smrtlink  test1  #将用户liyubing加入用户组smrtlink
# chmod  -R 775 /opt #smrtlink组需要权限在opt下创建文件夹


#以liyubing用户安装
# SMRT_ROOT=/data1/software/pacbio/smrtlink
# ./smrtlink_5.1.0.26412.run --rootdir $SMRT_ROOT --batch --ignore-system-check


##docker (ok)
docker pull jarvice/app-pacbio-smrtlink
docker run -it jarvice/app-pacbio-smrtlink /bin/bash

 #bash
SMRT_ROOT=/opt/pacbio/smrtlink
$SMRT_ROOT/admin/bin/services-start 
$SMRT_ROOT/admin/bin/services-status #(Optional)

