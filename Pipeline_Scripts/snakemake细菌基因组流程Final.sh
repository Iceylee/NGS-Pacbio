#测试目录： BYZH服务器
 /home/liyubing/analysis/7_snakemake/Bac_Genome/
#测试时间：2.6M基因组 10：18-10：52 

# 自己帐号测试
#（conda新建一个python3的环境，并安装snakemake）
conda create -n py3 python=3
source activate py3

conda install -c bioconda snakemake 
#第一次跑流程会先安装R依赖包

# liyubing帐号直接测试
ssh liyubing@192.168.31.156  
#password:liyubing
source activate py3

#文件说明
#rawData文件夹 基因组fasta（已组装） keg文件（kegg注释用）
#scripts 分析脚本
#rules 流程规则
#envs 脚本所需环境
#config.yaml 配置文件 
#snakefile-mine.py snakemake脚本
#color.txt circos作图着色用

#prop
1.将基因组fasta放入rawData文件夹并重命名为assembly.fasta （rawData中有的几个fasta文件是不同项目的细菌基因组，都可以用来测试。）
2.修改config.yaml配置文件。
fasta：基因组所在路径
name：基因组命名（一般以物种名拉丁文命名）
size2：基因组大小*2并取整（比如2.6M，这里写5）
contig：单个与否？TRUE则单个contig
db：需要注释的库去掉前面的井号
short_chr:circos图中每个contig的命名


#开始分析流程
sudo service docker start
# sudo usermod -a -G docker liyubing (帐号已加入组中用这个命令 安全性更高)
snakemake --snakefile snakefile-mine.py --use-conda

##报错
prokka报错原因：
tbl2asn过期。这个软件每隔一段时间会过期一次，如果遇到报错，下载最新的tbl2asn替换到prokka的binaries目录下即可。
可参考https://github.com/tseemann/prokka/issues/303
