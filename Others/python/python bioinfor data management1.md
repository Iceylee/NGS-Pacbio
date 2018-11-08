
1.给定FASTA格式的文件(test1.fa 和 test2.fa)，写一个程序 cat.py 读入文件，并输出到屏幕
* 用到的知识点
    * open(file)
    * for .. in loop
    * print
    * the amazng , or strip() function



```python
%%writefile script/cat.py

#作业1和2, cat.py
#包含了三种处理换行符的写法
for line in open("data/test1.fa"):
    print (line.strip())

# for line in open("data/test2.fa"):
#     print line,

# for line in open("data/test1.fq"):
#     line = line.strip()
#     print line
```

    Overwriting script/cat.py



```python
%run script/cat
```

    >NM_001011874 gene=Xkr4 CDS=151-2091
    gcggcggcgggcgagcgggcgctggagtaggagctggggagcggcgcggccggggaaggaagccagggcg
    >NM_001195662 gene=Rp1 CDS=55-909
    AGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCA
    >NM_0112835 gene=Rp15 CDS=128-6412
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCAC
    >NM_011283 gene=Rp1 CDS=128-6412
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCAC


二. 写程序 `splitName.py`, 读入test2.fa, 并取原始序列名字第一个空格前的名字为处理后的序列名字，输出到屏幕
    * 用到的知识点
        * split
        * 字符串的索引
    * 输出格式为：
```
>NM_001011874
gcggcggcgggcgagcgggcgctggagtaggagctg.......
```


```python
%%writefile script/splitName.py
for line in open("data/test2.fa"):
    if line[0] == '>':
        line = line.split(' ')[0]
    print (line.strip())
```

    Overwriting script/splitName.py



```python
%run script/splitName
```

    >NM_001011874
    gcggcggcgggcgagcgggcgctggagtaggagctggggagcggcgcggccggggaaggaagccagggcg
    aggcgaggaggtggcgggaggaggagacagcagggacaggTGTCAGATAAAGGAGTGCTCTCCTCCGCTG
    CCGAGGCATCATGGCCGCTAAGTCAGACGGGAGGCTGAAGATGAAGAAGAGCAGCGACGTGGCGTTCACC
    CCGCTGCAGAACTCGGACAATTCGGGCTCTGTGCAAGGACTGGCTCCAGGCTTGCCGTCGGGGTCCGGAG
    >NM_001195662
    AAGCTCAGCCTTTGCTCAGATTCTCCTCTTGATGAAACAAAGGGATTTCTGCACATGCTTGAGAAATTGC
    AGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCA
    AGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGT
    GGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGC
    TGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACA
    CAGCATCACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTG
    >NM_011283
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCAC
    ACAAACTATTTACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACAC
    CTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATAT
    CACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGG
    GTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCC
    TGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGA
    GGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAAGGCCCGCAGG
    CGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATA
    TGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA
    >NM_0112835
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCAC
    ACAAACTATTTACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACAC
    CTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATAT
    CACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGG
    GTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCC
    TGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGA
    GGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAAGGCCCGCAGG
    CGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATA
    TGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA



三. 写程序 `formatFasta.py`, 读入test2.fa，把每条FASTA序列连成一行然后输出
    * 用到的知识点
        * join
        * strip    
    * 输出格式为:
```
    >NM_001011874
    gcggcggcgggc......TCCGCTG......GCGTTCACC......CGGGGTCCGGAG
```


```python
%%writefile script/formatFasta.py

#读取一个fasta file，建立一个dictionary

#打开fasta文件
try:
    f = open("data/test2.fa")
except IOError:
    print ("File not exist!")
    
#建立dic
seqs = {}
for line in f:
    line = line.strip()#discard每行的 \n 符号
    if line[0] == '>': #如果是新的序列，该行为id
        words = line.split() #按空格分开
        name = words[0][1:] #第一个字段，去掉大于号
        seqs[name] = '' #字段为key，创建dic的条目
    else: #仍然在序列中
        seqs[name] = seqs[name] + line  #继续延长序列
f.close()

for name in seqs:
    print (name)
    print (seqs[name])

```

    Overwriting script/formatFasta.py



```python

```


```python
%run script/formatFasta
```

    NM_011283
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCACACAAACTATTTACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAAGGCCCGCAGGCGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATATGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA
    NM_001011874
    gcggcggcgggcgagcgggcgctggagtaggagctggggagcggcgcggccggggaaggaagccagggcgaggcgaggaggtggcgggaggaggagacagcagggacaggTGTCAGATAAAGGAGTGCTCTCCTCCGCTGCCGAGGCATCATGGCCGCTAAGTCAGACGGGAGGCTGAAGATGAAGAAGAGCAGCGACGTGGCGTTCACCCCGCTGCAGAACTCGGACAATTCGGGCTCTGTGCAAGGACTGGCTCCAGGCTTGCCGTCGGGGTCCGGAG
    NM_0112835
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCACACAAACTATTTACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAAGGCCCGCAGGCGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATATGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA
    NM_001195662
    AAGCTCAGCCTTTGCTCAGATTCTCCTCTTGATGAAACAAAGGGATTTCTGCACATGCTTGAGAAATTGCAGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTG


四. 写程序 `formatFasta-2.py`, 读入test2.fa，把每条FASTA序列分割成80个字母一行的序列
    * 用到的知识点
        * 字符串切片操作
        * range
    * 输出格式为
```    
    >NM_001011874
    gcggcggcgc.(60个字母).TCCGCTGACG #(每行80个字母)
    acgtgctacg.(60个字母).GCGTTCACCC
    ACGTACGATG(最后一行可不足80个字母)
```     


```python
%%writefile script/formatFasta-2.py
filename = "data/test2.fa" #将参数都提前写出来，方便统一修改
length = 80
seqList = {}

try:
    f = open(filename)
except IOError:
    print ("file not exists!")
    
for line in f:
    line = line.rstrip()
    if line[0] == ">":
        name = line.split(' ')[0]
        seqList[name] = ''
    else:
        seqList[name] += line

for name in seqList:
    print (name)
    seq = seqList[name]
    print (len(seq))
    for i in range(0,len(seq),length): #索引从0开始
        print (seq[i:i+length]) #每隔80切片 
```

    Overwriting script/formatFasta-2.py



```python
%run script/formatFasta-2.py
```

    >NM_011283
    630
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCACACAAACTATT
    TACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATG
    ATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAG
    TTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTC
    TGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATC
    ACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAA
    GGCCCGCAGGCGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATA
    TGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA
    >NM_0112835
    630
    AATAAATCCAAAGACATTTGTTTACGTGAAACAAGCAGGTTGCATATCCAGTGACGTTTATACAGACCACACAAACTATT
    TACTCTTTTCTTCGTAAGGAAAGGTTCAACTTCTGGTCTCACCCAAAATGAGTGACACACCTTCTACTAGTTTCTCCATG
    ATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATTCAAATATCACTCATCCTGTAGTGGCTAAACGCATCAG
    TTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTGGTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTC
    TGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGTAAGGAACATCAGCACGCCCCGTGGACGACACAGCATC
    ACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCTCCCACAATAAGAAGGTGCTGCCAGTTGACCTGGACAA
    GGCCCGCAGGCGCCCTCGGCCCTGGCTGAGTAGTCGCTCCATAAGCACGCATGTGCAGCTCTGTCCTGCAACTGCCAATA
    TGTCCACCATGGCACCTGGCATGCTCCGTGCCCCAAGGAGGCTCGTGGTCTTCCGGAATGGTGACCCGAA
    >NM_001195662
    420
    AAGCTCAGCCTTTGCTCAGATTCTCCTCTTGATGAAACAAAGGGATTTCTGCACATGCTTGAGAAATTGCAGGTCTCACC
    CAAAATGAGTGACACACCTTCTACTAGTTTCTCCATGATTCATCTGACTTCTGAAGGTCAAGTTCCTTCCCCTCGCCATT
    CAAATATCACTCATCCTGTAGTGGCTAAACGCATCAGTTTCTATAAGAGTGGAGACCCACAGTTTGGCGGCGTTCGGGTG
    GTGGTCAACCCTCGTTCCTTTAAGACTTTTGACGCTCTGCTGGACAGTTTATCCAGGAAGGTACCCCTGCCCTTTGGGGT
    AAGGAACATCAGCACGCCCCGTGGACGACACAGCATCACCAGGCTGGAGGAGCTAGAGGACGGCAAGTCTTATGTGTGCT
    CCCACAATAAGAAGGTGCTG
    >NM_001011874
    280
    gcggcggcgggcgagcgggcgctggagtaggagctggggagcggcgcggccggggaaggaagccagggcgaggcgaggag
    gtggcgggaggaggagacagcagggacaggTGTCAGATAAAGGAGTGCTCTCCTCCGCTGCCGAGGCATCATGGCCGCTA
    AGTCAGACGGGAGGCTGAAGATGAAGAAGAGCAGCGACGTGGCGTTCACCCCGCTGCAGAACTCGGACAATTCGGGCTCT
    GTGCAAGGACTGGCTCCAGGCTTGCCGTCGGGGTCCGGAG


五. 写程序 `sortFasta.py`, 读入test2.fa, 并取原始序列名字第一个空格前的名字为处理后的序列名字，排序后输出
    * 用到的知识点
        * sort
        * dict
        * aDict[key] = []
        * aDict[key].append(value)


```python
%%writefile script/sortFasta.py

filename = "data/test2.fa"
seqList = {}

with open(filename) as f:
    for line in f:
        line = line.strip()
        if line[0] == '>':
            name = line.split(' ')[0]
            seqList[name] = ''
        else:
            seqList[name] += line

#dictionary 按key排序方法。只能导出keys，sort之后，再打印
keys = seqList.keys()
keys.sort()

for key in keys:
    print (key)
    print (seqList[key])
    
            
```

    Overwriting script/sortFasta.py



```python
%run script/sortFasta.py
```


    ------------------------------------------------------------------------

    AttributeError                         Traceback (most recent call last)

    /Users/Icey/Downloads/NGS-Study/learn/ipython/script/sortFasta.py in <module>()
         14 #dictionary 按key排序方法。只能导出keys，sort之后，再打印
         15 keys = seqList.keys()
    ---> 16 keys.sort()
         17 
         18 for key in keys:


    AttributeError: 'dict_keys' object has no attribute 'sort'



```python

```


```python

```


```python

```


```python

```




    

        
5. 提取给定名字的序列
    * 写程序 `grepFasta.py`, 提取fasta.name中名字对应的test2.fa的序列，并输出到屏幕。
    * 写程序 `grepFastq.py`, 提取fastq.name中名字对应的test1.fq的序列，并输出到文件。
        * 用到的知识点
            * print >>fh, or fh.write()
            * 取模运算，4 % 2 == 0
    
6. 写程序 `screenResult.py`, 筛选test.expr中foldChange大于2的基因并且padj小于0.05的基，可以输出整行或只输出基因名字
    * 用到的知识点
        * 逻辑与操作符 and 
        * 文件中读取的内容都为字符串，需要用int转换为整数，float转换为浮点数

6. 写程序 `transferMultipleColumToMatrix.py` 将文件(multipleColExpr.txt)中基因在多个组织中的表达数据转换为矩阵形式
    * 用到的知识点
        * aDict['key'] = {}
        * aDict['key']['key2'] = value
        * if key not in aDict
        * aDict = {'ENSG00000000003': {"A-431": 21.3, "A-549", 32.5,...},"ENSG00000000003":{},}
    * 输入格式(只需要前3列就可以)
```
Gene    Sample  Value   Unit    Abundance
ENSG00000000003 A-431   21.3    FPKM    Medium
ENSG00000000003 A-549   32.5    FPKM    Medium
ENSG00000000003 AN3-CA  38.2    FPKM    Medium
ENSG00000000003 BEWO    31.4    FPKM    Medium
ENSG00000000003 CACO-2  63.9    FPKM    High
ENSG00000000005 A-431   0.0     FPKM    Not detected
ENSG00000000005 A-549   0.0     FPKM    Not detected
ENSG00000000005 AN3-CA  0.0     FPKM    Not detected
ENSG00000000005 BEWO    0.0     FPKM    Not detected
ENSG00000000005 CACO-2  0.0     FPKM    Not detected
```
    * 输出格式
```
Name	A-431	A-549	AN3-CA	BEWO	CACO-2
ENSG00000000460	25.2	14.2	10.6	24.4	14.2
ENSG00000000938	0.0	0.0	0.0	0.0	0.0
ENSG00000001084	19.1	155.1	24.4	12.6	23.5
ENSG00000000457	2.8	3.4	3.8	5.8	2.9
```

6. 写程序 `reverseComplementary.py`计算序列 `ACGTACGTACGTCACGTCAGCTAGAC`的反向互补序列
    * 用到的知识点
        * reverse
        * list(seq)
7. 写程序 `collapsemiRNAreads.py`转换smRNA-Seq的测序数据
    * 输入文件格式(mir.collapse, tab-分割的两列文件，第一列为序列，第二列为序列被测到的次数)
```
        ID_REF        VALUE
        ACTGCCCTAAGTGCTCCTTCTGGC        2
        ATAAGGTGCATCTAGTGCAGATA        25
        TGAGGTAGTAGTTTGTGCTGTTT        100
        TCCTACGAGTTGCATGGATTC        4
```
    * 输出文件格式 (mir.collapse.fa, 名字的前3个字母为样品的特异标示，中间的数字表示第几条序列，是序列名字的唯一标示，第三部分是x加每个reads被测到的次数。三部分用下划线连起来作为fasta序列的名字。)
```        
        >ESB_1_x2
        ACTGCCCTAAGTGCTCCTTCTGGC
        >ESB_2_x25
        ATAAGGTGCATCTAGTGCAGATA
        >ESB_3_x100
        TGAGGTAGTAGTTTGTGCTGTTT
        >ESB_4_x4
        TCCTACGAGTTGCATGGATTC
```
8. 简化的短序列匹配程序 (map.py) 把short.fa中的序列比对到ref.fa, 输出短序列匹配到ref.fa文件中哪些序列的哪些位置
    * 用到的知识点
        * find
    * 输出格式 (输出格式为bed格式，第一列为匹配到的染色体，第二列和第三列为匹配到染色体序列的起始终止位置（位置标记以0为起始，代表第一个位置；终止位置不包含在内，第一个例子中所示序列的位置是(199,208](前闭后开，实际是chr1染色体第199-206的序列，0起始). 第4列为短序列自身的序列.)。
    * 附加要求：可以只匹配到给定的模板链，也可以考虑匹配到模板链的互补链。这时第5列可以为短序列的名字，第六列为链的信息，匹配到模板链为'+'，匹配到互补链为'-'。注意匹配到互补链时起始位置也是从模板链的5'端算起的。
``` 
    chr1	199	208	TGGCGTTCA
    chr1	207	216	ACCCCGCTG
    chr2	63	70	AAATTGC
    chr3	0	7	AATAAAT
```

10. 备注：
    * 每个提到提到的“用到的知识点”为相对于前面的题目新增的知识点，请综合考虑。此外，对于不同的思路并不是所有提到的知识点都会用着，而且也可能会用到未提到的知识点。但是所有知识点都在前面的讲义部分有介绍。
    * 每个程序对于你身边会写的人来说都很简单，因此你一定要克制住，独立去把答案做出，多看错误提示，多比对程序输出结果和预期结果的差异。
    * 学习锻炼“读程序”，即对着文件模拟整个的读入、处理过程来发现可能的逻辑问题。
    * 程序运行没有错误不代表你写的程序完成了你的需求，你要去查验输出结果是不是你想要的。
11. 关于程序调试
    * 在初写程序时，可能会出现各种各样的错误，常见的有缩进不一致，变量名字拼写错误，丢失冒号，文件名未加引号等，这时要根据错误提示查看错误类型是什么，出错的是哪一行来定位错误。当然，有的时候报错的行自身不一定有错，可能是其前面或后面的行出现了错误。
    * **用脑袋运行程序**：当程序写作完成后，自己尝试对着数据文件，一行一行的执行程序，来看程序的运行是否与自己想干的活一致，有没有纰漏。
    * 当结果不符合预期时，要学会**使用print来查看每步的操作是否正确**，比如我读入了字典，我就打印下字典，看看读入的是不是我想要的，是否含有不该存在的字符；或者**在每个判断句、函数调入的情况下打印个字符，来跟踪程序的运行轨迹**。


```python
%%writefile script/formatFasta2.py

#读取一个fasta file，建立一个dictionary

#打开fasta文件
try:
    f = open("data/assembly.pep")
except IOError:
    print ("File not exist!")
    
#建立dic
seqs = {}
for line in f:
    line = line.strip()#discard每行的 \n 符号
    
    if line[0] == '>': #如果是新的序列，该行为id
        #words = line.split() #按空格分开
        name = line #第一个字段，去掉大于号
        seqs[name] = '' #字段为key，创建dic的条目
    else: #仍然在序列中
        seqs[name] = seqs[name] + line  #继续延长序列
f.close()

f1 = open("assembly.fasta","w")

for name in seqs:
    f1.writelines(name)
    f1.writelines("\n")
    f1.writelines (seqs[name])
    f1.writelines("\n")
    
f1.close()
  
```

    Overwriting script/formatFasta2.py



```python
%run script/formatFasta2.py
```


```python
!less assembly.fasta
```

    >unitig_0|quiver_3042 # 3437601 # 3438059 # 1 # ID=1_3042;partial=00;start_type=[m ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.686[m
    MNILGSHMSLDALHLKISSGMVVAARHWRRLCQGALTGYGISEACAVPLLMIVRLGDGVHQVAVAQAAGLESPSLVRLLD[m QLCKAGLVCRSEDPLDRRAKALSLTVEGRALAESIEGELVRLRREVLGGIDQADLDATLRVIQAFEAAGVMP*[m
    >unitig_0|quiver_4070 # 4582964 # 4584229 # 1 # ID=1_4070;partial=00;start_type=[m ATG;rbs_motif=AGGAG;rbs_spacer=5-10bp;gc_cont=0.607[m
    MGEKNIIVDRGIPPSGSSGGGGRGGSSTGGITIPLGVNGEPSHHALNVAALMNGTVLEAVLEGQGWPSPDAYYDLGVDMW[m GMLPYQIVEMRNELSDSYLRKERNLPASLNAELAAAEAAAGSTAALPDSKKAERSIGIIKSMMATRDQQIAFNRGRLATE[m QGGRFDDRSIKEVIDELRKLDDYDVPAALDVELSLYTAALALHVDLKAQEQLREKLDALEKARRDALEKESYKEAATYAS[m DIGKEIANRFGNQVAQAANDMQKGIAGKRIGSYQEALKAFEKLSQNPGLPLNAKDSVAVAQALEAPDKATLGDNMLRLGK[m AFGVTGGVIQAAGLVDSAVSGFKTGDWKPFLLELESAVVGKVAGSLAGAMVGIALGFLVSVPAGAVAGTVLAAVFIGAAS[m SYFDTEKVDQINQWVTGVVAP*[m
    >unitig_0|quiver_1132 # 1258627 # 1260837 # -1 # ID=1_1132;partial=00;start_type[m =ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.660[m
    MNATTRDNGGLRQRDLDALARAEHADPFAVLGPHDDGAGGLLVRAFLPNALNARLLARHDGQVLAEMVQGSVPGLFTAHL[m DQARAYLLQIGWAGGEQVTEDPYSFGPQLGDMDLYLFAEGNHRDLSGRFGAQPTVVDGVAGVCFSVWAPNARRVSVVGDF[m NNWDGRRHPMRLRHSAGVWELFVPRLGVGETYKFEVLGADGVLPLKADPLARATELPPSTASKVAGALAHDWRDGQWMEQ[m RAQRHAYSAPLSIYELHAGSWRCELDDAGEVARFYNWRELAERLVPYVQELGFTHIELMPIMEHPFGGSWGYQPLSLFAP[m TSRYGSAEDFAFFVDACHQGGIGVILDWVPAHFPTDEHGLARFDGTALYEYDNPLEGFHQDWNTLIYNLGRNEVRGFMLA[m SALHWLKHFHIDGLRVDAVASMLYRDYSRKAGEWVPNRHGGRENLEAIDFIRHLNGVAAHEAPGALIIAEESTAWPGVSQ[m PVQQGGLGFAYKWNMGWMHDTLHYIQNDPVHRTYHHNEMSFGLIYAYSEHFILPISHDEVVHGKHSLIDKMPGDRWQKFA[m NLRAYLAFMWTHPGKKLLFMGCEFGQWREWNHDQQLDWYLLQYSEHQGVQKLVADLNRLYRELPALHEQDCRAQGFQWLI[m GDDAHNSVYAWLRWSSQGEPLLVVANFTPVPREGYRIGVPFGERWQELLNSDAGLYAGSNVGNLGGVACEAIASHGQPLS[m [7massembly.fasta[m[K


```python

```
