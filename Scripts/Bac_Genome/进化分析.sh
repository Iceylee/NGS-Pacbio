#KA KS
#准备同源基因的pep和nuc文件。pep文件包含一对同源基因的蛋白序列。nuc包含该对的核苷酸序列。
#用clustalw2 比对这两条氨基酸序列 aln
clustalw2 two.pep

pal2nal.pl two.aln two.nuc -nogap >two.codon

#pal2nal.pl支持多种输出格式，如果下游要用PAML继续作分析，可以指定输出文件符合paml的输入格式：
pal2nal.pl two.aln two.nuc -output paml  -nogap  >  two.codon

#PAML需要3个输入文件：two.tree, two.codon, two.cnt (control file)
codeml two.cnt
#two.codeml


#kaks-calculator (如果序列完全一致，会一直执行不出结果)
pal2nal.pl ../two.aln ../two.nuc -nogap >two.codon
perl ~/bin/codon2kaks.pl two.codon >two.axt
KaKs_Calculator -i two.axt -o two.axt.kaks
KaKs_Calculator -i two.axt -o two.axt.kaks -m [LWL,YN,MYN]
#-m 参数可以指定各种计算ka ks的方法。
#-o 参数，指定输出文件名字。

#假基因
#test_parameter_file 修改脚本所在的路径
#fasty36路径：/home/liyubing/biosoft/fasta-36.3.8g/bin
#需要蛋白序列，核苷酸序列，gff文件(可以没有？)
 python ~/biosoft/PseudogenePipeline/_wrapper_scripts/CombinedPseudoWrapper.py test_parameter_file










