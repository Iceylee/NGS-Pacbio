
#islandpath install
sudo service docker start
docker pull brinkmanlab/islandpath:1.0.0


docker run \
	-it --rm \
	-v /home/liyubing/analysis/7_snakemake/Bac_Genome/:/data  brinkmanlab/islandpath:1.0.0 /apps/islandpath/example/NC_003210.gbk /data/NC_003210_gis.txt


#casFinder install

#ubuntu中装软件。从容器运行软件。
docker pull ubuntu

#apt-get问题 宿主机
nmcli dev show | grep 'IP4.DNS'
sudo vim /etc/docker/daemon.json
#加入下面信息
{
    "dns": ["192.168.31.1","8.8.8.8"]
}
#
#重启Docker服务，命令：
sudo service docker restart

#指定容器名字CasFinder
docker run -it --name CasFinder ubuntu
#再次进入这个容器
docker start CasFinder
docker attach CasFinder
#复制文件
docker cp -r /data1/software/CRISPRCasFinder CasFinder:/opt

#sudo command not found
su -
apt-get install sudo

#报错后缺少的一些东西安装
apt-get install apt-utils
apt-get install dialog 
apt-get install whiptail
apt-get install make
sudo apt-get install python-dev gcc

#
cpanm JSON::Parse --force

#test
CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def General -cas -i install_test/sequence.fasta -out Results_test_install –keep

#保存容器为新的镜像
docker ps -l #查看容器ID
docker commit e360ce2b52ee icey_casfinder #镜像名字要小写
docker images 


docker run -it --rm icey_casfinder /bin/bash
