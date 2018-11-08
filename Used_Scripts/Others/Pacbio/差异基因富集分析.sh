1.差异基因序列进行注释    
emapper.py -i DiffProtein_Fa/114_Diff.fa --output 114_maNOG -m diamond --cpu 16

2.获取其中的GO注释结果 
python /data1/script/GO_KEGG_Annotation/GetGOID_Annotation.py 114_maNOG.emapper.annotations  /data1/NCBI/gene2go  > 114_GO_annotation.txt

3.根据差异基因注释结果 以及 该物种的GO注释结果进行富集分析

python  /data1/script/GO_KEGG_Annotation/Pyher_GOEnrichment_Analysis.bk.py MF.txt all.MF.du.txt test > GO_FourValue.txt


#
cat MF.txt|cut -f 1|sort -k1,1 -u|wc -l
cat all.MF.du.txt|cut -f 1|sort -k1,1 -u|wc -l

