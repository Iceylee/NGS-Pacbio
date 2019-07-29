#无参
%s/\/data2\/ClientData\/2018_07\/ZhangYiWen\//\/data2\/ClientData\/2018_08\/DongWuZi\//g

python TransSeq_NoGenome.py TransSeq_NoGenome.conf 1>run.log 2>&1

tail -f run.log


#Result
mkdir Result && cd Result
mkdir 0.QC
mkdir 1.Trinity
mkdir 2.Mapping
mkdir 3.ExpressAnalysis
mkdir 4.DiffExprGene
mkdir 5.Annotation
mkdir 6.EnrichmentAnalysis
mkdir 7.TransDecoder
mkdir 8.CorrelationAnalysis

cd Result/
mv ../FastQC/* 0.QC

cd Result/1.Trinity
ln -s ../../1.Trinity/Trinity.fasta ../../1.Trinity/Trinity_N50.txt ./

cd ../2.Mapping
ln -s ../../2.Mapping/*bam_stat.out ../../2.Mapping/Statistic_Mapping.txt ./

cd ../3.ExpressAnalysis
ln -s ../../3.ExpressAnalysis/RSEM* ./

cd ../4.DiffExprGene
ln -s ../../4.DiffExprGene/DESeq2/3.DiffExprGene/* ./

cd ../5.Annotation
cp ../../5.Annotation/*.out ../../5.Annotation/*.txt ./ 

cd ../6.EnrichmentAnalysis
cp ../../6.EnrichmentAnalysis/* ./

cd ../7.TransDecoder
cp ../../7.TransDecoder/Trinity_CD-HIT_0.9.fa.transdecoder.* ./

#Annotation colnames
cat NR_Annotation.txt |sed "1 i Gene_ID\tNR_ID\tIdentity\tAlign_Length\tMismatch\tGap\tQ_start\tQ_end\tT_start\tT_end\tE_value\tScore\tFunction"  > temp
mv temp NR_Annotation.txt
cat KOG_Annotation.txt|sed "1 i Gene_ID\tKOG_ID\tIdentity\tAlign_Length\tMismatch\tGap\tQ_start\tQ_end\tT_start\tT_end\tE_value\tScore\tFunction"  > temp
mv temp KOG_Annotation.txt
cat SwissProt_Annotation.txt|sed "1 i Gene_ID\tSP_ID\tIdentity\tAlign_Length\tMismatch\tGap\tQ_start\tQ_end\tT_start\tT_end\tE_value\tScore\tFunction"  > temp
mv temp SwissProt_Annotation.txt

cat GO_Annotation.txt|sed "1 i Gene_ID\tGO_ID\tOntology\tTerm" > temp
mv temp GO_Annotation.txt

#anno 

awk '
BEGIN{FS=OFS="\t"}
NR==FNR{a[$1]=$13;next}
$0~/>/{split($0,x,">");
       split(x[2],y," ");
       p=y[1];print $0" "a[p]}
$0!~/>/{print $0}' Annotation_nr.out  Trinity_CD-HIT_0.9.fa >Trinity_CD-HIT_0.9_anno.fa



BamDir = OutPut + "1.Mapping/"
python /data1/script/CountTrinityBamState.py $BamDir
cd %s

for i in $(find . -name '*mapping.log')
do 
echo `dirname $i|awk -F'/' '{print $2}'`,`basename $i` >> mapping_stat_list.txt 
done 

for i in $(find . -name '*_stat.out')
do 
echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> bam_stat_list.txt 
done

python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > Statistic_Mapping.txt


#参考基因组
mkdir Genome
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.fa Genome/
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.gtf Genome/

#
cp /data1/script/Pipeline/RNA-Seq/* ./

#
mkdir R_input
mv colData.csv R_input/colData.csv.temp
mv sample_list.txt sample_list.txt.temp
mv Gene2Symbol.txt Gene2Symbol.txt.temp

ExtraExonLenFromGTF.py

#提取
awk '{                                
    for (i = 1; i <= NF; i++) {
        if ($i ~ /gene_id|gene_name/) {
            printf "%s ", $(i+1)
        }
    }
    print ""
}' *.gtf | sed -e 's/"//g' -e 's/;//g' -e 's/ /\t/' | sort -k1,1 | uniq > Gene2Symbol.txt

awk '{                                
    for (i = 1; i <= NF; i++) {
        if ($i ~ /gene_id|gene_biotype/) {
            printf "%s ", $(i+1)
        }
    }
    print ""
}' *.gtf | sed -e 's/"//g' -e 's/;//g' -e 's/ /\t/' | sort -k1,1 | uniq > Gene2Type.txt

cat Gene2Type.txt|grep "lincRNA"|cut -f 1 > LncRNA_ID.txt

%s/\/data3\/ClientData\/2018_07\//\/data2\/ClientData\/2018_08\//g

python /data1/script/ExtraExonLenFromGTF.py *.gtf > exonLen.txt

Rscript /data1/script/deseq2+GO+KEGG/Rpipe/0Get_OrgDb_KEGG_Name.R “Danio rerio”

#[1] "AH57981"
#   kegg_code scientific_name common_name
# 71       dre     Danio rerio   zebrafish

#snpeff
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.fa sequences.fa
cp /data1/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.gtf genes.gtf

java -jar snpEff.jar build -gtf22 -v Danio_rerio #gtf
#

#下载snpeff database
java -jar snpEff.jar databases > db.txt
java -jar snpEff.jar download Danio_rerio

cat genes.gtf.bk|grep -v -E "KN150086.1|KN150447.1|KN150272.1" >genes.gtf

#修改名字
for smp in WR181061S WR181062S
do
	mv ${smp}-1_clean_R1.fastq.gz ${smp}_1_clean_R1.fastq.gz
	mv ${smp}-1_clean_R2.fastq.gz ${smp}_1_clean_R2.fastq.gz
done


#seq lengths
/data1/software/trinityrnaseq-Trinity-v2.4.0/util/misc/fasta_seq_length.pl ../1.Trinity/CD-HIT/Trinity_CD-HIT_0.9.fa  > Trinity_CD-HIT_0.9.seqLengths

#anno stat
for file in Annotation_kog.out Annotation_nr.out Annotation_uniprot.out GO_Annotation.txt KEGG_Annotation.txt
do
	cat $file|sort -k1,1 -u|wc -l
done

33710
54528
41467
16724
25279
