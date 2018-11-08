#software

docker pull wwliao/cnvnator

docker run \
  --user root \
  -v /home/liyubing/analysis/test:/data \
  --rm \
  wwliao/cnvnator \
  cnvnator -root file.root -tree WR180001S_sorted.bam -unique

#bowtie
docker pull biodckrdev/bowtie:1.1.2
##在镜像里执行命令 数据输出到指定目录
docker run \
  --user root \
  -v /home/liyubing/analysis/test:/data \
  --rm \
  biodckrdev/bowtie:1.1.2 \
  bowtie-build /data/22.fa 22.fa


#常用
docker images  
docker ps #正在运行的容器
docker ps -a #所有已经创建的容器

#了解数据卷 volume
docker volume create hello
docker run --mount type=volume,source=hello,target=/data wwliao/cnvnator -unique -root out.root -tree WR180001S_sorted.bam
docker run -d --name nginx -p 80:80 -v /data/www/html:/usr/share/nginx/html nginx

#在容器中建立了test文件夹？
docker run -it --name container-test -h CONTAINER -v /test biodckrdev/bowtie:1.1.2 /bin/bash


#run images
docker run -it miarmaseq/centos /bin/bash --rm



#install software inside images？？

#data processing and then exit and then come back？？

#access data outside of container



#
docker pull busybox
docker run busybox ls #列出busybox镜像的所有内容
docker run busybox ls /bin #busybox的bin目录下的内容
uname -a #当前运行系统
docker run busybox uname -a








docker build -t friendlyhello .  # Create image using this directory's Dockerfile
docker run -p 4000:80 friendlyhello  # Run "friendlyname" mapping port 4000 to 80
docker run -d -p 4000:80 friendlyhello         # Same thing, but in detached mode
docker container ls                                # List all running containers
docker container ls -a             # List all containers, even those not running
docker container stop <hash>           # Gracefully stop the specified container
docker container kill <hash>         # Force shutdown of the specified container
docker container rm <hash>        # Remove specified container from this machine
docker container rm $(docker container ls -a -q)         # 删除所有容器
docker image ls -a                             # List all images on this machine
docker image rm <image id>            # Remove specified image from this machine
docker image rm $(docker image ls -a -q)   # Remove all images from this machine
docker login             # Log in this CLI session using your Docker credentials
docker tag <image> username/repository:tag  # Tag <image> for upload to registry
docker push username/repository:tag            # Upload tagged image to registry
docker run username/repository:tag                   # Run image from a registry




Volume只有在下列情况下才能被删除：
该容器是用docker rm －v命令来删除的（-v是必不可少的）。
docker run中使用了--rm参数







