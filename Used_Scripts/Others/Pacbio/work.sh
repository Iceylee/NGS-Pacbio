# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/675/GCF_000150675.1_ASM15067v2/GCF_000150675.1_ASM15067v2_genomic.fna.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/675/GCF_000150675.1_ASM15067v2/GCF_000150675.1_ASM15067v2_genomic.gff.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/400/815/GCA_000400815.2_VDAG_JR2v.4.0/GCA_000400815.2_VDAG_JR2v.4.0_genomic.fna.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/400/815/GCA_000400815.2_VDAG_JR2v.4.0/GCA_000400815.2_VDAG_JR2v.4.0_genomic.gff.gz
# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/400/815/GCA_000400815.2_VDAG_JR2v.4.0/GCA_000400815.2_VDAG_JR2v.4.0_genomic.gbff.gz

# wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-37/gff3/verticillium_dahliaejr2/Verticillium_dahliaejr2.VDAG_JR2v.4.0.37.gff3.gz

#http://fungi.ensembl.org/Verticillium_dahliae/Info/Index

wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-37/fasta/verticillium_dahliaejr2/cds/Verticillium_dahliaejr2.VDAG_JR2v.4.0.cds.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/fungi/release-37/fasta/verticillium_dahliae/cds/Verticillium_dahliae.ASM15067v2.cds.all.fa.gz

makeblastdb -in ASM15067v2.fa -dbtype nucl -out db/ASM15067v2

blastn -query VDAG_JR2.fa -db db/ASM15067v2 -out  results.txt -evalue 1e-5 -outfmt 6 -perc_identity 80 -max_target_seqs 1 -num_threads 8

#EGY14128	VDAG_05292
cat ASM15067v2.fa |grep '>'|awk 'BEGIN{FS=" ";OFS="\t"}{split($4,y,":");split($1,x,">");print x[2],y[2]}' > ASM15067v2.gene.tab
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$2;next}{print $1, a[$2]}' ASM15067v2.gene.tab results.txt > JR2toASM.tab
#去重1列 2列完全相同的重复
cat JR2toASM.tab|awk '!a[$1,$2]++' >JR2toASM.tab2

#重复的2列 提取
cat JR2toASM.tab2 |awk '{print $2}'|sort |uniq -d > duplicated.VDAG
#对应的第一列 
awk 'BEGIN{FS=OFS="\t"}FNR==NR{a[$1]=$1;next}($2 in a){print $1,$2}' duplicated.VDAG JR2toASM.tab2|sort -k1,1 -u | head