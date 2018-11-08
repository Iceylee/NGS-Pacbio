#### Variant filtering ####

variant_name=LM.samtools.vcf
filter_name=LM.fit1.vcf
filter_log_name=LM.flt1.log
GenomeFa=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa
GenomeDict=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.dict

#java -jar /data1/software/picard/picard.jar CreateSequenceDictionary R=$GenomeFa O=$GenomeDict #1min

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T VariantFiltration  -R $GenomeFa -V $variant_name -window 35 -cluster 3 -filterName FSS -filter 'FS > 30.0' -filterName QD -filter 'QD < 2.0' -o $filter_name 2>$filter_log_name #3h


 #### SNP/INDEL ####
SNP_name=LM.SNPs.vcf
INDEL_name=LM.INDELs.vcf
filter_name=LM.fit1.vcf
GenomeFa=/data1/GenomicDatabases/Maize/GCF_000005005.2_B73_RefGen_v4_genomic.fa

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType SNP -o $SNP_name # 1h

java -jar /data1/software/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar -T SelectVariants -R $GenomeFa -V $filter_name -selectType INDEL -o $INDEL_name #1h

 #### Annotation ####
filter_name=LM.fit1.vcf
Annota_Html=LM.SnpEff.html
Annota_vcf=LM.SnpEff.vcf
snpEff_SpeciesName=Maize
java -Xmx10G -jar /data1/software/snpEff/snpEff.jar eff -s $Annota_Html -c /data1/software/snpEff/snpEff.config -v -ud 500 $snpEff_SpeciesName $filter_name > $Annota_vcf