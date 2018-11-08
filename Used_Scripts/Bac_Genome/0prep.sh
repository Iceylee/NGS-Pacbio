#blast db

#1.db dir 加入环境变量bashrc
nrdir=~/database/nr/nr
trembldir=~/database/trembl/uniprot_trembl.fasta
pfamdir=~/database/pfam/Pfam-A-duplicate.fasta
spdir=~/database/swissport/uniprot_sprot.fasta
cogdir=~/database/cog/prot2003-2014.fa
vfdbadir=~/database/VFdb/VFDB_setA_pro.fas
vfdbbdir=~/database/VFdb/VFDB_setB_pro.fas
ardbdir=~/database/ardbAnno/blastdb/resisGenes.pfasta
cazydir=~/database/CAzy/CAZydb
phidir=~/database/phi/phi

#2.index
makeblastdb -in CAZyDB.07202017.fa -dbtype prot -out CAZydb 
#加-parse_seqids报错 
# BLAST Database creation error: Error: Duplicate seq_ids are found:
# LCL|AHJ76089.2|GH73


#3.anno
#从database的fa文件获取所有序列信息列，得到文件.anno。存入database/anno
cat database/cog/prot2003-2014.fa | grep ">" > ~/analysis/pro1_pseudomonas_putida/blast/cog.anno
#nr.anno得到后，需要再次处理，只保留 \x01\分隔符的前面部分
cat nr.anno|awk '{FS='\x01';print $1}' > nr.anno2
mv nr.anno2 nr.anno


#4.analysis/database

#4.1.gaf处理
#下载goa.gaf , 根据第2和5列，去重复
grep -v "\!" goa_uniprot_all.gaf |awk '!a[$2,$5]++' > goa_uniprot_all_duplicates_moved.gaf
#提取第二列spID和第5列GO，
#将相同的spID合并。第一列相同的合并第二列，分号间隔
cat goa_uniprot_all_duplicates_moved.gaf|awk '{FS=OFS="\t"}{print $2,$5}'>sp_go.gaf
awk -f merge.awk sp_go.gaf > sp_go_merge.gaf

#4.2.pfam
#提取pfamACC(只取点号前）和term并去重复(database)
cat pdb_pfam_mapping.txt |awk 'BEGIN{OFS="\t"}{split($5,x,".");print x[1],$7}'|sort -k1,1 -u>pfam_term_dup_moved.tab
#4.3.cog
#将逗号转成tab分隔
cat cog2003-2014.csv |awk 'BEGIN{FS=",";OFS="\t"}{$1=$1}1' > cog2003-2014.csv.tab




echo 'BEGIN { OFS="\t" }
{
    curr = $1
    if (curr == prev) {
        rec = rec ";" $2
    }
    else {
        if (rec) print rec
        rec = $0
    }
    prev = curr
}
END { if (rec) print rec }' > merge.awk