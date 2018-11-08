NarrowPeak = OutPut + "2.CallPeak/" + Contrast + "_peaks.broadPeak"
            Peak_bed = OutPut + "4.GetPeakGene/" + Contrast + "_peaks.bed"
            Peak_Annotation = OutPut + "3.Annotation/" + Contrast + "_peaks_annootation.txt"
            Sample_Genesat3KborlessfromPeaks = OutPut + "4.GetPeakGene/" + Contrast + "_Genesat3KborlessfromPeaks.txt"
            Sample_closestGene = OutPut + "4.GetPeakGene/" + Contrast + "_closesstGeneAt3KborLess.bed"

            Command10_1 = "annotatePeaks.pl %s %s -gtf %s > %s" % (NarrowPeak, GG
enomebed, GenomeGTF, Peak_Annotation)
            print "--** annotation peak **--Command10_1: \n%s\n\n" % Command10_1
            os.system(Command10_1)

            Command10_2 = "less -S %s | cut -f1-4,9 > %s" % (NarrowPeak, Peak_bee
d)
            Command10_3 = "intersectBed -wa -a %s -b %s > %s" % (HumanGenesPlusMM
inus3kb, Peak_bed, Sample_Genesat3KborlessfromPeaks)
            print "--** annotation peak **--Command10_2: \n%s\n\n%s\n\n" % (Commm
and10_2, Command10_3)