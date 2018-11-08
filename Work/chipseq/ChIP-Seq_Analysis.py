#!/usr/bin/env python

### Author: zss
### Date: 2018.7.3
## For: ChIP-Seq with Genome ###

import ConfigParser
import sys
import os
import commands

if len(sys.argv) != 2:
    print "\n\t\tUsage: python %s ChIP-Seq_Analysis.conf\n\n" % sys.argv[0]
    sys.exit(1)

ConfgFile = sys.argv[1]

Config = ConfigParser.ConfigParser()
Config.read(ConfgFile)

### Get IP Info
SamplesID = Config.get('IP', 'SamplesID')
#Control = Config.get('IP', 'Control')
TreatInfo = eval(Config.get('IP', 'TreatInfo'))

### Get DB Info
GenomeFa = Config.get("DB", "GenomeFa")
GenomeGTF = Config.get("DB", "GenomeGTF")
GeneNames = Config.get("DB", "GeneNames")
Genomebed = Config.get('DB', "Genomebed")
GenomeName = Config.get('DB', 'GenomeName')
HumanGenesPlusMinus3kb = Config.get('DB', 'HumanGenesPlusMinus3kb')

### Get Data Info
InPut = Config.get("Data", "InPut")
OutPut = Config.get("Data", "OutPut")

### SNP
GATK = Config.get("SNP", "GATK")
Human = Config.get("SNP", "Human")

species_name = Config.get("SNP", "species_name")

#####################################################################################################################
print "\n################################################ Info ###############################################################\n"
print "\t\tGenome File: %s"  % GenomeFa
print "\t\tData Dir: %s"  % InPut
print "\t\tOutPut Dir: %s"  % OutPut 
print "\n#####################################################################################################################\n"
#####################################################################################################################


### File + "_R1.fastq.gz" or  + "_R2.fastq.gz"
def GetCleanData(InPut):
    List = []
    FileList = os.listdir(InPut)
    for i in FileList:
        if i.endswith("_R1.fastq.gz"):
            SeqName = i.split("_R1.fastq.gz")[0]
            if SeqName not in List:
                List.append(SeqName)
    print "\t\t\033[1;31;40m *** Sample: %s \033[0m \n\n" % ('\t'.join(List))
    return List

def MkdirOutPut(OutPut):
    for i in [ "0.QC", "1.Mapping", "2.CallPeak", "3.Annotation/distance_to_TSS", "3.Annotation/region_classify", "3.Annotation/coverage", "4.GetPeakGene", "5.FindMotif", "6.DiffPeak", "7.SNP", "5.SNP_Indel", "8.OtherAnalysis/PearsonCorr", "8.OtherAnalysis/FingerPrint", "8.OtherAnalysis/TSS_Heatmap", "R_input"]:
        if not os.path.exists(OutPut + i):
            os.makedirs(OutPut + i)

def Mapping(InPut, OutPut, FileList, GenomeFa, GenomeGTF):
    genome_index = GenomeFa.replace(".fa", ".4.bt2")
    if not os.path.exists(genome_index):
        Command0 = "bowtie2-build --threads 15 %s %s" % (GenomeFa, GenomeFa.replace(".fa", ""))
        print "\033[1;31;40m --** bowtie2-build **--Command0: === %s ===   \033[0m \n\n" % (Command0)
        os.system(Command0)

    for i in FileList:
        Fq_R1 = InPut + i + "_R1.fastq.gz"
        Fq_R2 = InPut + i + "_R2.fastq.gz"

        #print "Fastq File: \n%s\n%s\n" % (Fq_R1, Fq_R2)
        MappingOut = OutPut + "1.Mapping" + "/" + i + ".bam"
        MappingOut_sorted_bam = OutPut + "1.Mapping" + "/" + i + "_sorted.bam"
        Bowtie2_log = OutPut + "1.Mapping" + "/" + i + "_bowtie2.log"
        HtseqOut   = OutPut + "1.Mapping" + "/" + i + "_CountNum.txt"

        ## 1. Mapping - bowtie
        Command1 = "bowtie2 -p 15 -x %s -1 %s -2 %s 2>%s | samtools view -bS - > %s" % (GenomeFa.replace(".fa", ""), Fq_R1, Fq_R2, Bowtie2_log, MappingOut)
        print "\033[1;31;40m --** Mapping **--Command1: === %s ===   \033[0m \n%s\n\n" % (i, Command1)
        os.system(Command1)

        ## 2. Samtools sort and deduplication and index 
        Command2 = "samtools sort -o %s %s -@ 10 " % (MappingOut_sorted_bam, MappingOut)
        print "\033[1;35;40m --** Sam2SortedBam **--Command2: \033[0m \n%s\n\n" % Command2
        os.system(Command2)

        Command3 = "samtools view -q 20 %s > %s " % (MappingOut_sorted_bam, MappingOut)
        print "\033[1;35;40m --** Deduplication **--Command3: \033[0m \n%s\n\n" % Command3
        os.system(Command3)

        Command4 = "samtools index %s " % (MappingOut_sorted_bam)
        print "\033[1;35;40m --** Index **--Command4: \033[0m \n%s\n\n" % Command4
        os.system(Command4)


### Sample Correlation Plot
def GetReadsCountMatrix_Plot(OutPut, SamplesID):
    MappingOut = OutPut + "1.Mapping/"
    ReadCount_npz = OutPut + "1.Mapping/readCounts.npz"
    ReadCount_tab = OutPut + "1.Mapping/readCounts.tab"

    ## 3. Get reads count matrix 
    Command5 = "multiBamSummary bins --bamfiles %s*_sorted.bam --minMappingQuality 30  --labels %s -out %s --outRawCounts %s" % (MappingOut, SamplesID, ReadCount_npz, ReadCount_tab)
    print "\033[1;35;40m --** Count Matrix **--Command5: \033[0m \n%s\n\n" % Command5
    os.system(Command5)

    ## 4. plotCorrelation
    scatterplot_png = OutPut + "8.OtherAnalysis/PearsonCorr/scatterplot_PearsonCorr_bigwigScores.png"
    PearsonCorr_tab = OutPut + "1.Mapping/PearsonCorr_bigwigScores.tab"
    Command6 = "plotCorrelation -in %s --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Samples'  --whatToPlot scatterplot -o %s --outFileCorMatrix %s" % (ReadCount_npz, scatterplot_png, PearsonCorr_tab)
    print "\033[1;35;40m --** plotCorrelation: scatterplot  **--Command6: \033[0m \n%s\n\n" % Command6
    os.system(Command6)

    heatmap_png = OutPut + "8.OtherAnalysis/PearsonCorr/heatmap_SpearmanCorr_readCounts.png"
    SpearmanCorr_tab = OutPut + "1.Mapping/SpearmanCorr_readCounts.tab"
    Command7 = "plotCorrelation -in %s --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Samples' --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o %s --outFileCorMatrix %s" % (ReadCount_npz, heatmap_png, SpearmanCorr_tab)
    print "\033[1;35;40m --** plotCorrelation: heatmap  **--Command7: \033[0m \n%s\n\n" % Command7
    os.system(Command7)


### FingerPrint
def FingerPrint(OutPut, SamplesID):
    Bam = OutPut + "1.Mapping/*_sorted.bam"
    FingerPrint_png = OutPut + "8.OtherAnalysis/FingerPrint/fingerprints.png"
    FingerPrint_tab = OutPut + "8.OtherAnalysis/FingerPrint/fingerprints.tab"

    ## 5. fingerprint 
    Command8 = "plotFingerprint -b %s --labels %s --minMappingQuality 30 --skipZeros -T 'Fingerprints of different samples' --plotFile %s --outRawCounts %s" % (Bam, SamplesID, FingerPrint_png, FingerPrint_tab)
    print "\033[1;35;40m --** fingerprint: **--Command8: \033[0m \n%s\n\n" % Command8
    os.system(Command8)


### Get Bam States
def GetBamState(OutPut):
    BamDir = OutPut + "1.Mapping/"
    mapping_stat = OutPut + "1.Mapping/mapping_stat_list.txt"
    bam_stat = OutPut + "1.Mapping/bam_stat_list.txt"
    Statistic_out = OutPut + "1.Mapping/Statistic_Mapping.txt"

    Command3_1 = "python /data1/script/CountBamState.py %s\n\n" % (BamDir)
    print "--** CountBamState **--Command3: \n%s\n\n" % Command3_1
    Command5_1 = "cd %s; for i in $(find . -name '*_bowtie2.log');do echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> %s ;done"  % (BamDir, mapping_stat)
    Command5_2 = "cd %s; for i in $(find . -name '*_stat.out');do echo `basename $i|awk -F'_' '{print $1}'`,`basename $i` >> %s;done"  % (BamDir, bam_stat)
    Command5_3 = "cd %s; python /data1/script/GetBamStat.py bam_stat_list.txt mapping_stat_list.txt > %s" % (BamDir, Statistic_out)
    os.system(Command3_1)
    os.system(Command5_1)
    os.system(Command5_2)
    os.system(Command5_3)


### Plot TSS Heatmap
def Plot_TSS_Heatmap(OutPut, Genomebed):
    BamDir = OutPut + "1.Mapping/"
    for i in os.listdir(BamDir):
        if i.endswith("_clean_sorted.bam"):
            bigWig = OutPut + "1.Mapping/" + i.replace("_clean_sorted.bam", "_clean_sorted.bigWig")
            Command6 = "cd %s; bamCoverage --bam %s -o %s" % (BamDir, i, bigWig)
            print "--** bam to bigWig **--Command6: \n%s\n\n" % Command6
            os.system(Command6)

    List = []
    Matrix_Enhancers = OutPut + "8.OtherAnalysis/TSS_Heatmap/matrix_Enhancers.tab.gz"
    Matrix_Enhancers_Dir = OutPut + "8.OtherAnalysis/TSS_Heatmap"
    for i in os.listdir(BamDir):
        if i.endswith("_clean_sorted.bigWig"):
            List.append(i)
    print "-----------+++\n%s\n%s----------------+++" % (List, ' '.join(List))
    Command7 = "cd %s; computeMatrix reference-point -S %s -R %s --referencePoint TSS -b 3000 -a 3000 -out %s" % (BamDir, ' '.join(List), Genomebed, Matrix_Enhancers)
    print "--** bigWig to matrix **--Command7: \n%s\n\n" % Command7
    os.system(Command7)
  
    Command8 = "cd %s; bash /data1/script/Pipeline/ChIP-Seq/TSS_Heatmap_Plot.sh %s" % (Matrix_Enhancers_Dir, Matrix_Enhancers)
    print "--**  matrix to TSS Heatmap **--Command8: \n%s\n\n" % Command8
    os.system(Command8)
    


### Call Peak and Annotation && ChIPseeker 
def CallPeak(OutPut, TreatInfo, GenomeName, Genomebed, GenomeGTF, HumanGenesPlusMinus3kb):
#def CallPeak(OutPut, Control, GenomeName, Genomebed, GenomeGTF, HumanGenesPlusMinus3kb):
    
    Peak_dir = OutPut + "2.CallPeak/"
    ### TreatInfo = {"ChLib28": ["ChLib26", "ChLib27"]}  # CK:[Treat1, Treat2]
    for i in TreatInfo:
        for j in TreatInfo[i]:
            Input = OutPut + "1.Mapping/" + i + "_clean_sorted.bam"
            IP = OutPut + "1.Mapping/" + j + "_clean_sorted.bam"
            Contrast = j + "_vs_" + i
            Contrast_log =  Peak_dir + Contrast + ".masc2.log"

            ## Macs2 callpeak
            Command9 = "macs2 callpeak -t %s -c %s -m 5 50 -p 1e-5 -f BAM -g %s -n %s --outdir %s -B --broad 2>%s" % (IP, Input, GenomeName, Contrast, Peak_dir, Contrast_log)
            print "--** callpeak **--Command9: \n%s\n\n" % Command9
            os.system(Command9)

            ## Get Peaks Gene and Annotation
            NarrowPeak = OutPut + "2.CallPeak/" + Contrast + "_peaks.broadPeak"
            Peak_bed = OutPut + "4.GetPeakGene/" + Contrast + "_peaks.bed"
            Peak_Annotation = OutPut + "3.Annotation/" + Contrast + "_peaks_annotation.txt"
            Sample_Genesat3KborlessfromPeaks = OutPut + "4.GetPeakGene/" + Contrast + "_Genesat3KborlessfromPeaks.txt"
            Sample_closestGene = OutPut + "4.GetPeakGene/" + Contrast + "_closestGeneAt3KborLess.bed"

            Command10_1 = "annotatePeaks.pl %s %s -gtf %s > %s" % (NarrowPeak, Genomebed, GenomeGTF, Peak_Annotation)
            print "--** annotation peak **--Command10_1: \n%s\n\n" % Command10_1
            os.system(Command10_1)

            Command10_2 = "less -S %s | cut -f1-4,9 > %s" % (NarrowPeak, Peak_bed)
            Command10_3 = "intersectBed -wa -a %s -b %s > %s" % (HumanGenesPlusMinus3kb, Peak_bed, Sample_Genesat3KborlessfromPeaks)
            print "--** annotation peak **--Command10_2: \n%s\n\n%s\n\n" % (Command10_2, Command10_3)
            os.system(Command10_2)
            os.system(Command10_3)

            ## Get Peaks Annotation: coverage, region_classify, distance_to_TSS
            NarrowPeak_Chr = OutPut + "2.CallPeak/" + Contrast + "_peaks_Chr.broadPeak"
            Peak_Annotation_Dir = OutPut + "3.Annotation/"
            Command10_5 = "python /data1/script/MacsFormatChange.py addChr %s > %s" % (NarrowPeak, NarrowPeak_Chr)
            Command10_6 = "Rscript /data1/script/PeakAnnotation.R %s %s %s" % (NarrowPeak_Chr, Contrast, Peak_Annotation_Dir)
            print "--** ChIPseeker Annotation **--Command10_5/6: \n%s\n\n%s\n\n" % (Command10_5, Command10_6)
            os.system(Command10_5)
            os.system(Command10_6)


### Get Motif
def GetMotif(OutPut, TreatInfo, GenomeFa):
    for i in TreatInfo:
        for j in TreatInfo[i]:
            Contrast = j + "_vs_" + i
            NarrowPeak = OutPut + "2.CallPeak/" + Contrast + "_peaks.broadPeak"
            NarrowPeak_fa = OutPut + "5.FindMotif/" + Contrast + "_peaks.fa"

            Command11 = "bedtools getfasta -fi %s -bed %s -name > %s" % (GenomeFa, NarrowPeak, NarrowPeak_fa)
            print "--** Get Peak Fa **--Command11: \n%s\n\n" % Command11
            os.system(Command11)

            NarrowPeak_homer_bed = OutPut + "5.FindMotif/" + Contrast + "_peaks_homer.bed"
            MotifDir = OutPut + "5.FindMotif/"
            Command12_1 = "NarrowPeak2HomerBed.sh %s > %s" % (NarrowPeak, NarrowPeak_homer_bed)
            Command12_2 = "findMotifsGenome.pl %s %s %s -len 8,10,12" % (NarrowPeak_homer_bed, GenomeFa, MotifDir)
            print "--** Get Motif **--Command12_1: \n%s\n\n" % Command12_1
            print "--** Get Motif **--Command12_2: \n%s\n\n" % Command12_2
            os.system(Command12_1)
            os.system(Command12_2)


### DiffBind Peak  and Diff Gene $|Fold| >= 2


### Get Diff Peak Motif


if __name__ == '__main__':
    print "\t\t\033[1;32;40m ### Step1: mkdir OutPut Dir ###\033[0m\n\n"
    #MkdirOutPut(OutPut)   

    print "\t\t\033[1;32;40m ### Step2: get fastq List ### \033[0m \n\n"
    FileList = GetCleanData(InPut)

    print "\t\t\033[1;32;40m ### Step3: Mapping & sort bam ###\033[0m \n\n"
    #Mapping(InPut, OutPut, FileList, GenomeFa, GenomeGTF)

    print "\t\t\033[1;32;40m ### Step4: State bam ###\033[0m \n\n"
    #GetBamState(OutPut)

    print "\t\t\033[1;32;40m ### Step4: Get Reads Count Matrix ###\033[0m \n\n"
    #GetReadsCountMatrix_Plot(OutPut, SamplesID)

    print "\t\t\033[1;32;40m ### Step5: FingerPrint ###\033[0m \n\n"
    #FingerPrint(OutPut, SamplesID)

    print "\t\t\033[1;32;40m ### Step6: Call Peak and Annotation ###\033[0m \n\n"
    CallPeak(OutPut, TreatInfo, GenomeName, Genomebed, GenomeGTF, HumanGenesPlusMinus3kb)

    print "\t\t\033[1;32;40m ### Step7: Get Motif ###\033[0m \n\n"
    GetMotif(OutPut, TreatInfo, GenomeFa)

    print "\t\t\033[1;32;40m ### Step8: Get TSS Heatmap ###\033[0m \n\n"
    Plot_TSS_Heatmap(OutPut, Genomebed)

    print "\t\t\033[1;32;40m ### Step9: Get SNP ###\033[0m \n\n"
    if GATK == "TRUE":
	from SNP_Indels import GetSNP_Indels
        #GetSNP_Indels(GenomeFa, OutPut, species_name, Human)
