

#TOFU
seqkit stat TOFU/collapse/tofu.collapsed.rep.fa #58,097

#matchanno
sed '1 d' matchAnnot_result.txt |wc -l #58072 稍微的差别：避免报错而去掉的少量输入序列
sed '1 d' matchAnnot_result.txt |grep "no_genes"|wc -l #840 没有找到gene
# gene数目
cat temp5|sort -k2,2 -u|wc -l #9529
# tr数目
cat temp5|sort -k3,3 -u|wc -l #12636

#基因结构优化
#总调整区域数（exon）
sed '1 d' extend_result.out|wc -l #33463
#reads
sed '1 d' extend_result.out|sort -k8,8 -u|wc -l #26062
# gene
cat exon.temp4|sort -k2,2 -u|wc -l #7045
# tr
cat exon.temp4|sort -k3,3 -u|wc -l #9010

#lncRNA
wc -l CPAT.list
wc -l CNCI.list
wc -l PLEK.list
#取交集
wc -l CNCI_PLEK_CPAT.list
#去掉ORF
seqkit stat lncRNA.final.fasta

#blast
wc -l nr.anno.out
wc -l sp.anno.out
wc -l trembl.anno.out




