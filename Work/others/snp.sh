/data1/software/annovar/humandb/hg38_avsnp150.txt


#dbSNP vcf 
ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-All.vcf.gz
#other version
ftp://ftp.ncbi.nih.gov/snp/organisms/
#human_9606_b151_GRCh38p7/	


#only annotate ID
java -jar /data1/software/snpEff/SnpSift.jar annotate -id /data1/software/annovar/humandb/00-All.vcf ChLib2.final.vcf > out.vcf

java -jar /data1/software/snpEff/SnpSift.jar annotate -id /data/dbsnp/All_20180418.vcf ChLib2.final.vcf >out.vcf #8.2 10am

#annotate ID & all INFO fields
java -jar SnpSift.jar annotate db/dbSnp/dbSnp137.20120616.vcf test.chr22.vcf 
