import sys,os
os.chdir("/Users/Icey/work/2018/0615snakemake")
input_list=[]
for i in units.itertuples():
    input_dir="htseq/%s_%s_CountNum.txt" % (i.sample, i.unit)
    input_list.append(input_dir)
    


#####count-matrix.py###

import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in input_list]

for t, (sample, unit) in zip(counts, units.index):
    t.columns = ["%s_%s" % (sample,unit)] #不懂框框是做什么 但里面放的是string


matrix = pd.concat(counts, axis=1) #横向合并。按匹配行名
matrix.index.name = "gene"

print(matrix)
matrix.to_csv("htseq_out.tab", sep="\t")


##star
import pandas as pd

counts = [pd.read_table(f, index_col=0, usecols=[0, 1], header=None, skiprows=4)
          for f in input_list]

#list，包含6个pandas df

for t, (sample, unit) in zip(counts, units.index):
    t.columns = [sample]

#t与counts对应；（sample，unit）与index中的sample，unit分别对应
#zip 把item与list中每一个对应。循环
matrix = pd.concat(counts, axis=1) #横向合并6个df
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1).sum() #按列名，相同的求和。
matrix.to_csv("all_counts.csv", sep="\t") #导出csv

总结：
1. 用于将多个样本经star输出的raw count合并成count matrix供后续deseq2分析
2. 每个tab只提取第一列gene名字和第二列（后面还有四列不知用途）。然后列合并。
3. 同一个样本的三个重复计算sum。最终的matrix只有两列count值。