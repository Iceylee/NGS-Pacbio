#!/bin/bash

# Choose Tophat OR Hisat
### mkdir 0.QC 1.Mapping 2.GenesExpress 3.DiffExprGene 4.GO_KEGG 5.SNP_Indel 6.l
ncRNA_Expr 7.Correlation_analysis

if [ $# != 4 ];then
        echo
        echo "          USAGE: $0 tophat/hisat DataDir DataGenomeFa OutputPath"
        echo "          e.g.: $0 hisat /data/ /data/hg38.genome OutputPath  (Out
putPath do not have '/'!)"
        echo
        exit 1;
fi

Path=$2
Genome=$3
OutputPath=$4

if [[ -d ${OutputPath} ]];then
        mkdir -p ${OutputPath}
fi
cd ${Path}
function GetFileList(){
        if [ $# != 1 ];then
                echo "pay attention the UASGE"
        fi
        Sample=()
        num=0
        for i in $(ls $1)
        do
                if [[ $i =~ 001.fq.gz ]] || [[ $i =~ 001.fastq.gz ]] ;then
                name=`echo $i|awk -F'_R' '{print $1}'`
                #name=`echo $i|awk -F'_' '{print $1"_"$2"_"$3"_"$4"_"$5}'`
                if ! [[ "${Sample[@]}" =~ $name ]];then
                        echo $name
                        Sample[num]=$name
                        num+=1
                fi
        fi
        done
        return ${Sample}
}

Sample=$(GetFileList $Path)
echo ${Sample}

### Choose gtf/gff/gff3
if [[ -f ${Genome}.gtf ]];then
        GenomeGTF=${Genome}.gtf
elif [[ -f ${Genome}.gff ]];then
        GenomeGTF=${Genome}.gff
elif [[ -f ${Genome}.gff3 ]];then
        GenomeGTF=${Genome}.gff3
fi
###


----------------------"

for i in ${Sample[@]}
do
        dirname=`echo $i|awk -F'_' '{print $1"_"$2}'`
        echo $dirname

        if [[ 'tophat' == $1 ]]
        then
                echo "Step1. Tophat Analysis! $i"
_R1_001.fastq.gz  ${i}_R2_001.fastq.gz

                echo "Step2. HTSeq-count Analysis! $i"
                if [[ -f ${GenomeGTF} ]]
                then
tOut/accepted_hits.bam ${GenomeGTF} -q > ${OutputPath}/${i}_count_number.txt &
                else
                        echo "Missed GTF/GFF file!"
                        continue
                fi

{dirname}_CuffOut ${dirname}_TophatOut/accepted_hits.bam

        elif [[ 'hisat' == $1 ]]
        then
                echo "Step1. Hisat Analysis! $i"
-2 ${i}_R2_001.fastq.gz -S ${OutputPath}/$i.sam

                echo "Step2. HTSeq-count Analysis! $i"
> ${OutputPath}/${i}_count_number.txt &

.sam

nd -o ${OutputPath}/${dirname}_CuffOut ${OutputPath}/${i}.bam &

gtf -l ${i} ${OutputPath}/${i}.bam
        else
ataGenomeFa "
        fi
done


--------------------"

echo "CountBamState!"
python /data/script/CountBamState.py ${OutputPath}



echo "Step3. RPKM Analysis! $i"
eneExonLen.txt

if [[ -f sample_list.txt ]];then
        mv sample_list.txt sample_list_bk.txt
fi

for i in ${Sample[@]}
do
    echo "$i    ${i}_count_number.txt" >> sample_list.txt
done

tPath}/sample_list.txt  > ${OutputPath}/AllSamplesRPKMValue.txt


----------"

echo "Step4. CountMatrix Analysis! $i"
python /data/script/CountMatrix.py sample_list.txt > CountMatrix4DESeq.csv

echo "------------------------------DONE-------------------------------------"