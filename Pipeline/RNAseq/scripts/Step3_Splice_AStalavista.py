#!/usr/bin/env python

# date: 2017/11/24 v1.0   # date: 2018/2/9 v2.0
# author: zss

### Input All transcripts.gtf list file
### Output Five type Count  "*_Alternative_summary.txt"

import sys,commands,os


def CountAstala(GTF_astalavista):
    Sample = "CK"
    Dict = {}
    Type1 = 0
    Type2 = 0
    Type3 = 0
    Type4 = 0
    Type5 = 0
    with open(GTF_astalavista, 'r') as F:
        for line in F:
            info = line.strip().split(';')
            structure = info[3].split('"')[1]
            #print structure
            if "0,1-2^" in structure:
                #print "0,1-2^: %s"  % structure
                Type1 += 1
            elif "1^,2^" in structure:
                Type2 += 1
            elif "1-,2-" in structure:
                Type3 += 1
            elif "0,1^2-" in structure:
                #print "0,1^2-: %s"  % structure
                Type4 += 1
            else:
                #print "other: %s" % structure
                Type5 += 1

    #print "Type1: %s, Type2: %s, Type3: %s, Type4: %s, Type5: %s" % (Type1, Type2, Type3, Type4, Type5)
    OutPutInfo1 = "Type1-Cassette exon, Type2-Alternative 5' splice site, Type3-Alternative 3' splice site, Type4-Retained intron, Type5-Mutually exclusive exon"
    OutPutInfo2 =  str(Type1) + ',' + str(Type2) + ',' + str(Type3) + ',' + str(Type4) + ',' + str(Type5)
    Dict[Sample]  = {'Type1': Type1, 'Type2': Type2, 'Type3': Type3, 'Type4':Type4, 'Type5':Type5}
    return Dict



if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "\tUsage: python %s Sample_transcript.txt\n" % sys.argv[0]
        sys.argv[1]

    File = sys.argv[1]

   ### Count Astalavista type
    print "Sample, Type1-Cassette exon, Type2-Alternative 5' splice site, Type3-Alternative 3' splice site, Type4-Retained intron, Type5-Mutually exclusive exon"
    Dict = CountAstala(File)
    sample = Dict.keys()[0]
    print "%s, %d, %d, %d, %d, %d" % (sample, Dict[sample]['Type1'], Dict[sample]['Type2'], Dict[sample]['Type3'], Dict[sample]['Type4'], Dict[sample]['Type5'])

