nohup /data/software/diamond-linux64/diamond makedb --in Danio_rerio.GRCz10.fa -d Danio_rerio.GRCz10.fa -p 22 &


nohup /data/software/diamond-linux64/diamond makedb --in Danio_rerio.GRCz10.cdna.all.fa -d Danio_rerio.GRCz10.cdna.all.fa -p 20 &

nohup /data/software/diamond-linux64/diamond makedb --in Danio_rerio.GRCz10.pep.all.fa -d Danio_rerio.GRCz10.pep.all.fa -p 20 &

diamond=/data/software/diamond-linux64/diamond
dbdir=/data/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.fa.dmnd
fasta=cDNA.unique.fa
dbout=DR.fa.out
nohup $diamond blastx -d ${dbdir} -q $fasta -o $dbout -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &


diamond=/data/software/diamond-linux64/diamond
dbdir=/data/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.cdna.all.fa.dmnd
fasta=cDNA.unique.fa
dbout=DR.cdna.fa.out
nohup $diamond blastn -d ${dbdir} -q $fasta -o $dbout -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle &

########
diamond=/data/software/diamond-linux64/diamond
dbdir=/data/GenomicDatabases/Danio_rerio/Danio_rerio.GRCz10.pep.all.fa.dmnd
fasta=cDNA.unique.fa
dbout=DR.pep.fa.out
nohup $diamond blastx -d ${dbdir} -q $fasta -o $dbout -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle &

cat DR.pep.fa.out|sort -k1,1 -u > DR.pep.fa.sort.out



nohup blastn -d ${dbdir} -q $fasta -o $dbout -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle &

