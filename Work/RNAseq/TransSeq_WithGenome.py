#!/usr/bin/env python

### Author: zss
### Date: 2018.3.30
## For: RNA-Seq with Genome ###

import ConfigParser
import sys
import os
#import rpy2.robjects as robjects
from subprocess import Popen
import commands


if len(sys.argv) != 2:
    print "\n\t\tUsage: python %s TransSeq_WithGenome.conf\n" % sys.argv[0]
    sys.exit(1)

ConfgFile = sys.argv[1]

Config = ConfigParser.ConfigParser()
Config.read(ConfgFile)

### Get DB Info
GenomeFa = Config.get("DB", "GenomeFa")
GenomeGTF = Config.get("DB", "GenomeGTF")
GeneExonLen = Config.get("DB", "GeneExonLen")
GeneNames = Config.get("DB", "GeneNames")

### Get Data Info
DataDir = Config.get("Data", "Dir")
OutPut = Config.get("Data", "OutPut")

### Get SampleList Info
SampleListFile = Config.get("SampleList", "SampleListFile")
#print "SampleListFile: %s" % SampleListFile

### Get LncRNA Info
LncRNAID_File = Config.get("LncRNA", "LncRNAID_File")

### Get RepNum Info
#RepNum = Config.get("Rep", "RepNum")

### Get R Info
CountTable = Config.get("R", "CountTable")
ColDataFile = Config.get("R", "ColDataFile")
BaseGroup = Config.get("R", "BaseGroup")

### Module
Module = Config.get("Module", "Module")
Species_Num = Config.get("Module", "Species_Num")
Species_Abbr = Config.get("Module", "Species_Abbr")
GO_Gene_Type = Config.get("Module", "GO_Gene_Type")
KEGG_Gene_Type = Config.get("Module", "KEGG_Gene_Type")

### SNP
GATK = Config.get("SNP", "GATK")

#####################################################################################################################
print "\n################################################ Info ###############################################################\n"
print "\t\tGenome File: %s"  % GenomeFa
print "\t\tData Dir: %s"  % DataDir
print "\t\tOutPut Dir: %s"  % OutPut 
print "\n#####################################################################################################################\n"
#####################################################################################################################


### File + "_R1.fastq.gz" or  + "_R2.fastq.gz"
def GetCleanData(DataDir):
    List = []
    FileList = os.listdir(DataDir)
    for i in FileList:
        if i.endswith("_R1.fastq.gz"):
            SeqName = i.split("_R1.fastq.gz")[0]
            if SeqName not in List:
                List.append(SeqName)
    print "\t\t\033[1;31;40m *** Sample: %s \033[0m \n\n" % ('\t'.join(List))
    return List

def MkdirOutPut(OutPut):
    for i in [ "0.QC", "1.Mapping", "2.GenesExpress", "3.DiffExprGene", "4.GO_KEGG", "5.SNP_Indel", "6.lncRNA_Expr", "7.Correlation_analysis", "R_input"]:
        if not os.path.exists(OutPut + i):
            os.makedirs(OutPut + i)

def MappingCount(DataDir, OutPut, FileList, GenomeFa, GenomeGTF):
    for i in FileList:
        Fq_R1 = DataDir + i + "_R1.fastq.gz"
        Fq_R2 = DataDir + i + "_R2.fastq.gz"

        #print "Fastq File: \n%s\n%s\n" % (Fq_R1, Fq_R2)
        MappingOut = OutPut + "1.Mapping" + "/" + i + ".sam"
        MappingOut_sorted_bam = OutPut + "1.Mapping" + "/" + i + "_sorted.bam"
        Hisat2_log = OutPut + "1.Mapping" + "/" + i + "_hisat2.log"
        HtseqOut   = OutPut + "1.Mapping" + "/" + i + "_CountNum.txt"

        Command1 = "hisat2 -p 15 -x %s -1 %s -2 %s -S %s > %s 2>&1\n\n" % (GenomeFa.replace(".fa", ""), Fq_R1, Fq_R2, MappingOut, Hisat2_log)
        print "\033[1;31;40m --** Mapping **--Command1: === %s ===   \033[0m \n%s" % (i, Command1)
        os.system(Command1)
        #if State1 != 0:
        #    print "\n\n!!! Error!!! Hisat2 Failed!!\n"
        #    sys.exit(2)
            
        Command2 = "htseq-count -f sam %s %s -q > %s\n\n" % (MappingOut, GenomeGTF, HtseqOut)
        print "\033[1;34;40m --** HtSeq-count **--Command2: \033[0m \n%s" % Command2
        os.system(Command2)

        Command3 = "samtools view -u %s | samtools sort -@ 15 - > %s\n\n" % (MappingOut, MappingOut_sorted_bam)
        print "\033[1;35;40m --** Sam2SortedBam **--Command3: \033[0m \n%s" % Command3
        os.system(Command3)

        #if State3 == 0:
        #    os.system("rm -rf MappingOut")


def GetBamState(OutPut):
#def GetBamState(OutPut, RSeQC_bam_stat_List, Mapping_List):
    BamDir = OutPut + "1.Mapping/"
    Command4_1 = "python /data1/script/CountBamState.py %s\n\n" % (BamDir)
    print "--** CountBamState **--Command4: \n%s" % Command4_1
    Command5_1 = "cd %s; for i in $(find . -name '*_hisat2.log');do echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> mapping_stat_list.txt ;done"  % BamDir
    Command5_2 = "for i in $(find . -name '*_stat.out');do echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> bam_stat_list.txt ;done" 
    Command5_3 = "python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > %s/1.Mapping/Statistic_Mapping.txt" % OutPut
    os.system(Command4_1)
    os.system(Command5_1)
    os.system(Command5_2)
    os.system(Command5_3)

def RPKM_Analysis(GenomeGTF, GeneExonLen, SampleListFile, OutPut, GeneNames):
    RPKM_OutFile = OutPut + "2.GenesExpress/AllSamplesRPKMValue.txt"
    Matrix_OutFile = OutPut + "R_input/CountMatrix4DESeq.csv"

    Command5 = "/data1/script/AnalysisRPKM.py %s %s > %s\n\n" % (GeneExonLen, SampleListFile, RPKM_OutFile)
    print "--** Count FPKM **--Command5: \n%s" % Command5
    os.system(Command5)

 
# 如果conf文件里定义了GeneNames，就输出带symbol的基因名称，否则，只输出geneid号
    if os.path.isfile(GeneNames):
        Command6 = "python /data1/script/CountMatrix.py %s %s > %s\n\n" % (SampleListFile, GeneNames, Matrix_OutFile)
    else:
        Command6 = "python /data1/script/CountMatrix.py %s > %s\n\n" % (SampleListFile, Matrix_OutFile)
    print "--** GetCountMatrix **--Command6: \n%s" % Command6
    os.system(Command6)


### Get LncRNA Expression and Count
def GetLncRNAExp(LncRNAID_File, OutPut, BaseGroup):
    AllSamplesRPKMValue = OutPut + "2.GenesExpress/AllSamplesRPKMValue.txt"
    AllSamplesCount = OutPut + "R_input/CountMatrix4DESeq.csv"
    LncRNAExp = OutPut + "6.lncRNA_Expr/LncRNA_RPKMValue.txt"

    LncRNACount = OutPut + "6.lncRNA_Expr/CountMatrix4DESeq.csv"
    ColDataFile = OutPut + "R_input/colData.csv"
    LncRNA_OutDir = OutPut + "6.lncRNA_Expr/"

    Command8_1 = "python /data1/script/GetlncRNARPKM.py RPKM %s %s > %s\n\n" % (LncRNAID_File, AllSamplesRPKMValue, LncRNAExp)
    Command8_2 = "python /data1/script/GetlncRNARPKM.py DESeq %s %s > %s\n\n" % (LncRNAID_File, AllSamplesCount, LncRNACount)
    Command8_3 = "Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R %s %s %s %s \n\n" % (LncRNACount, ColDataFile, BaseGroup, LncRNA_OutDir)
    print "--** GetlncRNARPKM **--Command8_1: \n%s" % Command8_1
    print "--** GetlncRNARPKM **--Command8_2: \n%s" % Command8_2
    print "--** GetlncRNARPKM **--Command8_3: \n%s" % Command8_3
    os.system(Command8_1)
    os.system(Command8_2)
    os.system(Command8_3)
 

### DESeq2 
def Rscript(CountTable, ColDataFile, BaseGroup, OutPut, Species_Num, Species_Abbr, GO_Gene_Type, KEGG_Gene_Type, Module):

    if os.path.isfile(CountTable):
        print "\t\t\033[1;32;40m ### Step8: R DESeq2 & HeatMap ###\033[0m \n\n"
        print "%s %s %s %s" % (CountTable, ColDataFile, BaseGroup, OutPut)
        Command10 = "Rscript /data1/script/deseq2+GO+KEGG/Rpipe/run.R %s %s %s %s"  % (CountTable, ColDataFile, BaseGroup, OutPut)
        print "--** Rscript **--Command10: \n%s" % Command10
        #os.system(Command10)
    
    if Module == "TRUE":
        print "\t\t\033[1;32;40m ### Step9: R  GO KEGG ###\033[0m \n\n"
        Command11 = "Rscript /data1/script/deseq2+GO+KEGG/Rpipe/3GO-KEGG.R %s %s %s %s %s"  % (Species_Num, Species_Abbr, GO_Gene_Type, KEGG_Gene_Type, OutPut)
        print "--** Rscript **--Command10: \n%s" % Command11
        os.system(Command11)


####   print "\033[1;31;40m --** Mapping **--Command1: === %s ===   \033[0m \n%s" % (i, Command1)

if __name__ == '__main__':
    print "\t\t\033[1;32;40m ### Step1: mkdir OutPut Dir ###\033[0m\n\n"
    MkdirOutPut(OutPut)   

    print "\t\t\033[1;32;40m ### Step2: get fastq List ### \033[0m \n\n"
    FileList = GetCleanData(DataDir)

    print "\t\t\033[1;32;40m ### Step3: Mapping & SeqCount ###\033[0m \n\n"
    MappingCount(DataDir, OutPut, FileList, GenomeFa, GenomeGTF)

    print "\t\t\033[1;32;40m ### Step4: State bam ###\033[0m \n\n"
    GetBamState(OutPut)

    print "\t\t\033[1;32;40m ### Step5: RPKM Analysis ###\033[0m \n\n"

    if os.path.isfile(SampleListFile):
        print "see"
        RPKM_Analysis(GenomeGTF, GeneExonLen, SampleListFile, OutPut, GeneNames)
    else:
        print "\n\t\t Please Make Sure Your SampleListFile! \n"
        sys.exit(2)

    if os.path.isfile(LncRNAID_File):
        print "\t\t\033[1;32;40m ### Step6: LncRNA RPKM ###\033[0m \n\n"
        GetLncRNAExp(LncRNAID_File, OutPut, BaseGroup)
    else:
        print "\n\t\t This Species GTF has no lncRNA~ Skip LncRNA Exp Analysis\n"

    print "\t\t\033[1;32;40m ### Step7: GATK/Samtools SNP ###\033[0m \n\n"
    if GATK:
        from SNP_Indels import GetSNP_Indels
        GetSNP_Indels(GenomeFa, OutPut)
    #else:
    #    from SNP_Indels_samtools import GetSNP_Indels
        #GetSNP_Indels(GenomeFa, OutPut, RepNum)
        

    # Rscript Analysis DESeq2 
    # Rscript GO KEGG
    Rscript(CountTable, ColDataFile, BaseGroup, OutPut, Species_Num, Species_Abbr, GO_Gene_Type, KEGG_Gene_Type, Module)
    
