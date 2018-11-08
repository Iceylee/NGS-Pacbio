/data1/script/Step1_FastQC.sh

/data1/script/Step1_FastQC_analysis.sh  
#质控分析脚本  /bin/bash /root/script/Step1_FastQC_analysis.sh  /data/ClientData/WuZhiGuo_MianHua/



#质控

#所有文件名需要更改为*_R1.fastq.gz   *_R2.fastq.gz
for file in `ls *_1.fq.gz`
do
	filename=${file/_1.fq.gz/}
	echo "sudo mv ${filename}_1.fq.gz ${filename}_R1.fastq.gz"
	echo "sudo mv ${filename}_2.fq.gz ${filename}_R2.fastq.gz"
done

#数据目录/Databackup/2018_06/ChenJinMing  
#结果目录/data2/ClientData/2018_06/ChenJinMing    
bash /data1/script/Step1_FastQC_analysis_PE.sh /Databackup/2018_06/ChenJinMing  /data2/ClientData/2018_06/ChenJinMing

#比对+差异分析
/data1/script/Pipeline/TransSeq_WithGenome.py   
#数据目录 /data2/ClientData/2018_06/JiangYang2/clean_data


index=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91
outdir=/home/liyubing/analysis/9_JiangYang
datadir=/data2/ClientData/2018_06/JiangYang2/clean_data/
gtf=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf
cd $datadir
for file in `ls *_R1.fastq.gz`
do
	filename=${file/_R1.fastq.gz/}
	samfile=$outdir/1.Mapping/$filename.sam
	hisat2_log=$outdir/1.Mapping/${filename}_hisat2.log
	echo "hisat2 -p 15 -x $index -1 $datadir/${filename}_R1.fastq.gz -2 $datadir/${filename}_R2.fastq.gz -S $samfile > ${hisat2_log} 2>&1"
	echo "$filename hisat2 over"

	
done > $outdir/cmd.list
nohup bash cmd.hisat2.list &



outdir=/home/liyubing/analysis/9_JiangYang/
datadir=/data2/ClientData/2018_06/JiangYang2/clean_data/
gtf=/data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf

SampleListFile=/data2/ClientData/2018_06/JiangYang2/sample_list.txt
Matrix_OutFile=$outdir/R_input/CountMatrix4DESeq1.csv
Matrix_Names_OutFile=$outdir/R_input/CountMatrix4DESeq.csv
GeneNames=/data1/GenomicDatabases/Human/Ensembl/Human_GeneID2EntrezID.csv





Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_input/CountMatrix4DESeq1.csv R_input/colData.csv  EV ./
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_12_45/CountMatrix4DESeq.csv R_12_45/colData.csv  EV R_12_45/
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_12_456/CountMatrix4DESeq.csv R_12_456/colData.csv  EV R_12_456/

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R R_12_3456/CountMatrix4DESeq.csv R_12_3456/colData.csv  EV R_12_3456/

/home/liyubing/analysis/9_JiangYang/R_12_456
#去除531


#R_12_456
cd R_12_456
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R CountMatrix4DESeq.csv colData.csv  EV ./
Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R AH57973 hsa ENSEMBL ENTREZID ./


datadir=/data2/ClientData/2018_06/JiangYang2/1.Mapping/
stringdir=/data1/software/stringtie-1.3.3b.Linux_x86_64/
#10min /per
for i in $(ls *_sorted.bam); do
        
        foo=${i/_clean_sorted.bam/}
        
        echo "$stringdir/stringtie $datadir/$i -p 10 -o stringtie_out/${foo}.gtf -G /data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf -l $foo"
done > /home/liyubing/analysis/9_JiangYang/1.Mapping/stringtie.cmd

cd /home/liyubing/analysis/9_JiangYang/1.Mapping

stringdir=/data1/software/stringtie-1.3.3b.Linux_x86_64/
cd stringtie_out
#1min 
nohup $stringdir/stringtie --merge -G /data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf -p 8 -o merged_EV.gtf mergelist_EV.txt &
nohup $stringdir/stringtie --merge -G /data1/GenomicDatabases/Human/Ensembl/Homo_sapiens.GRCh38.91.gtf -p 8 -o merged_NFAT1.gtf mergelist_NFAT1.txt &

#2min
nohup /data1/software/astalavista-4.0/bin/astalavista -t asta --threads 10 -i merged_EV.gtf &
nohup /data1/software/astalavista-4.0/bin/astalavista -t asta --threads 10 -i merged_NFAT1.gtf &

gunzip merged_EV_sorted.gtf_astalavista.gtf.gz
gunzip merged_NFAT1_sorted.gtf_astalavista.gtf.gz

python Step3_Splice_AStalavista.py merged_EV_sorted.gtf_astalavista.gtf >EV_count.txt
python Step3_Splice_AStalavista.py merged_NFAT1_sorted.gtf_astalavista.gtf >NFAT1_count.txt


