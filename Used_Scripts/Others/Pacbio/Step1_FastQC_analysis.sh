#!/bin/bash

Path=$1
cd $Path

Sample=('')
num=0
for i in $(ls $Path)
do
	if [[ $i =~ fastq.gz ]];then
		name=`echo $i|awk -F'_R' '{print $1}'`
		if ! [[ "${Sample[@]}" =~ $name ]];then
			Sample[num]=$name
			num+=1
		fi
	fi
done

for i in ${Sample[@]}
do 
	echo $i
	#dirname=`echo $i|awk -F'_' '{print $2}'`
	echo "--------------Step1:  原始数据质控分析"
	mkdir FastQC_raw
	/opt/FastQC/fastqc -o FastQC_raw -t 10 ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz
	echo "--------------Step2:  原始数据过滤"
	mkdir clean_data
	java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar  PE -threads 15  ${i}_R1_001.fastq.gz ${i}_R2_001.fastq.gz clean_data/paired_${i}_R1_001.fastq.gz clean_data/unpaired_${i}_R1_001.fastq.gz clean_data/paired_${i}_R2_001.fastq.gz clean_data/unpaired_${i}_R2_001.fastq.gz  ILLUMINACLIP:/opt/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 

	echo "--------------Step3:  过滤后数据质控分析"
	mkdir FastQC_clean
	/opt/FastQC/fastqc -o FastQC_raw -t 5 clean_data/paired_${i}_R1_001.fq.gz clean_data/paired_${i}_R2_001.fq.gz
done

#### 分析原始数据的Q30
cd FastQC_raw
python /root/script/Q30_analy.py .

#### 分析过滤后数据的Q30
cd FastQC_clean
python /root/script/Q30_analy.py .
