###常用

#传输文件的命令
scp nohup.out byzh2:/data/ 

#换权限（root 到 liyubing）
#更改目录所属用户及用户组  -R 是递归的意思
sudo chown liyubing:bioinfo -R R_input/
sudo chmod g+w -R Dongwuzi/

#查看文件  
ls -lrt #以修改时间逆序显示  
ls -alt #修改时间正序显示  
less -SN

#服务器内存  
free -g
#磁盘空间  
df -h
#文件按从大到小排列
du --max-depth=1 -h

#查找文件    
ls *.png #查找当前目录下所有 png 文件，ls **/*.png递归查找。  
find / -name YAML.pm -type f 2> >(grep -v 'Permission denied' >&2)
find ./ -name ".*" -print #查找目录（多层）隐藏文件



#装perl包  
perl -MCPAN -e shell #（进入cpan终端安装）
>install Data::Dumper
#或者conda装 conda install perl-data-dumper



###awk sed等文档操作
#粘贴列
paste -d , date1.csv date2.csv 
 
#在文件第7行后插入一行。本文件修改
sed -i '7ichromosomes={params.short_chr}' ./plot/ticks.conf

#vim 替换多个空字符\s  
%s/\s*/

#显示各种特殊字符
cat -A file 
vim         set list

#数字排序 倒序
sort -t $'\t' -k5 -nr
#科学记数法排序 正序
sort -t $'\t' -k5 -g

#vim 中替换
:%s/data/outData/g
全局




###
#列出 除去某个文件的其他所有
list -I
大i（ignore）

#xargs
ls -I pro1|xargs -I {} mv {} pro1/

#
screen -R xxx    
screen -r xxx  
#返回
ctrl+A+D     
#退出
ctrl+D

screen -ls

tail -f file
#（实时显示）


#计算器  
`bc -l`
#解压  
`gzip -d xxx.gz`  
`tar -xf xxx.tar`
#下载  
`curl -O 地址`





grep -E -o

aa=`command`

#bash条件判断。注意空格，中括号，还有是等号而不是双等号
if [ $a = 5 ]; then
command1
else
command2
fi

#去掉第一个字符 大于号 
cut -c 2- #从第二个字符开始
cat assembly.fasta|grep ">"|cut -c 2-

#重命名
按照file指定新名字和旧名字，将文件名中相应字符进行替换（rename 注意centos版本和用法）。
filename中tab分隔。分别赋给了old和new两个变量
#rename
while read -r old_name new_name; do
    rename $old_name $new_name *$old_name*.txt
done < filename

#rename _h _half *.png    centos7
#rename 's/$old_name/$new_name/' *$old_name*.txt

