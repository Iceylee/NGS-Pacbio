source activate py3
pip install snakemake

snakemake -np

snakemake

#RNA-seq workflow

####1.raw data
cd snakemake           # go to your home directory (please adapt)
mkdir -p projects/t-cell_analysis/samples/raw # create a directory containing projects
cd projects/t-cell_analysis/samples/raw       # change directory

# Download the file of interest (here using a loop)
# Note that it can be interesting to store the STDERR of wget.
for i in control_R1 control_R2 treated_R1 treated_R2
do 
	echo "wget http://pedagogix-tagc.univ-mrs.fr/courses/data/ngs/abd/${i}.fq.gz"
done

###2.bowtie2 index
cd snakemake # please adapt
mkdir -p projects/indexes/bowtie2/mm10 # Create a directory to store the index
cd projects/indexes/bowtie2/mm10       # go to projects/indexes/bowtie2/mm10
# Download chr19 sequence (mm10 version)
wget --no-clobber http://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz 2> chr19.fa_wget.log
gunzip chr19.fa.gz # uncompress the file
# to get help about bowtie-build type:
# bowtie2-build -h
bowtie2-build chr19.fa chr19 &> chr19.fa_bowtie2-build.log
ls -rtlh

###3.download GTF
cd snakemake/projects/ # adpat to your needs
mkdir -p annotations/gtf
cd annotations/gtf
#提取以19开头的，并将19换成chr19
wget ftp://ftp.ensembl.org/pub/release-83/gtf/mus_musculus/Mus_musculus.GRCm38.83.chr.gtf.gz | gunzip -c | grep "^19" | sed 's/^19/chr19/' > GRCm38.83.chr19.gtf

###4.创建Snakefile
cd projects/t-cell_analysis
touch Snakefile

# Adapt the path to your needs
cd snakemake/projects/t-cell_analysis
snakemake --snakefile snakefile.py

#安装需要py3的环境。运行不需要

/data1/software/miniconda2/envs/py3/bin/snakemake

#http://pedagogix-tagc.univ-mrs.fr/courses/ABD/practical/snakemake/projects/t-cell_analysis/

##5.画图 
snakemake --dag|dot -Tsvg > dag.svg
snakemake --dag 2> /dev/null | dot -T png > ../../img/workflow_bysample.png
snakemake --use-conda --snakefile snakefile.py



##python内引用R script需要安装rpy2
conda install rpy2

snakemake --cores 10



rnaseq-star-deseq2

snakemake -n
###star index
cd /home/liyubing/analysis/7_snakemake/rna-seq-star-deseq2/data/
STAR --runThreadN 15 --runMode genomeGenerate --genomeDir index/ --genomeFastaFiles refs/22.fa  --sjdbGTFfile refs/22.gtf --sjdbOverhang 100

snakemake --report report.html
