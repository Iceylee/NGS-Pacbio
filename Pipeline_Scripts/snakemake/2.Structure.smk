
rule repeatmasker:
  input: config["fasta"]
  output: "outData/repeatmasker/RepeatMasker.txt"
  params: threads=config["threads"],outdir="outData/repeatmasker/"
      
  shell: """
         mkdir -p {params.outdir}

         /data1/software/RepeatMasker/RepeatMasker -parallel {params.threads} -html -gff -dir {params.outdir} {input}
		
		name=`ls {params.outdir}/*.out`
		bash scripts/rm.sh $name {output}


  """

rule trf:
  input: config["fasta"]
  output: "outData/trf/Trf.txt"
  params: outdir="outData/trf/",workdir=config["workdir"]
      
  shell: """
  		 mkdir -p {params.outdir} && cd {params.outdir}

         trf {params.workdir}/{input} 2 7 7 80 10 50 500 -f -d -m
         cp *.dat out.dat

         awk '{{print $0 "Sequence: "> "file" NR}}' RS='Sequence: ' out.dat

         for file in `ls file*`
		 do
		 	grep -E '^[0-9]|unitig' $file|awk 'BEGIN{{FS=" ";OFS="\t";}}NR==1{{name=$1}}{{$1=$1}}NR>1{{print name"\t"$0}}' > $file.out
		 done
		 cat *.out >> all.out
		
		 cat all.out|cut -f 1,2,3,4,5,6,7,8,9,14,15,16|sed "1 i Sequence\tStart\tEnd\tPeriodSize\tCopyNumber\tConsensusSize\tPercentMatches\tPercentIndels\tScore\tEntropy\tConsensusSequences\tRepeatSequences" > {params.workdir}/{output}

  """



rule sRNA:
  input: config["fasta"]
  output: out="outData/sRNA/cmscan.out",tb="outData/sRNA/cmscan.tblout"
  params: threads=config["threads"],clanin=config["Rfam"]["clanin"],cm=config["Rfam"]["cm"],Z=config["size2"]
      
  shell: """

         cmscan -Z {params.Z} --cut_ga --rfam --nohmmonly --tblout {output.tb} --fmt 2 --cpu {params.threads} --clanin {params.clanin} {params.cm} {input} > {output.out}


  """

rule sRNA_stat:
  input: "outData/sRNA/cmscan.tblout"
  output: stats="outData/sRNA/ncRNA.stats",length="outData/sRNA/ncRNA.length",gff="outData/sRNA/sRNA.gff"
  params: anno=config["Rfam"]["anno"],outdir="outData/sRNA/"
      
  shell: """

         bash scripts/sRNA.sh {input} {params.anno} {output.stats} {output.length} {output.gff} {params.outdir}


  """


rule tRNA:
  input: config["fasta"]
  output: out="outData/tRNA/tRNA.out",ss="outData/tRNA/tRNA.ss",stats="outData/tRNA/tRNA.stats",gff="outData/tRNA/tRNA.gff"   
  shell: """

         tRNAscan-SE --bact -o {output.out} -f {output.ss} -m {output.stats} {input}
		 
    		 cat {output.out}|awk '
    		 BEGIN{{FS=OFS="\t"}}
    		 NR>3{{print $1,"tRNAscan-SE","tRNA",$3,$4,$9,"+",$5,$6}}'|sed "1 i #Sequence\tSource\tFeature\tStart\tEnd\tScore\tStrand\tType\tCodon" > {output.gff}

  """



rule rRNA:
  input: config["fasta"]
  output: fna="outData/rRNA/rRNA.fna",report="outData/rRNA/rRNA.hmmreport",xml="outData/rRNA/rRNA.xml",out="outData/rRNA/out.gff",gff="outData/rRNA/rRNA.gff"
  shell: """

         perl /data1/software/rnammer/rnammer -S bac -m lsu,ssu,tsu -multi -f {output.fna} -h {output.report} -xml {output.xml} -gff {output.out} {input}

         grep -v "#" {output.out}|awk 'BEGIN{{OFS="\t"}}{{print $1,"rRNAmmer",$3,$4,$5,$6,$7,$8,$9}}' |sed "1 i #Sequence\tSource\tFeature\tStart\tEnd\tScore\tStrand\tFrame\tAttribute" > {output.gff}		 

  """



rule GI:
  input: config["fasta"]
  output: out="outData/island/all.GIs.out",island="outData/island/islands.txt",GeneIsland="outData/island/GeneIslands.txt"
  params: outdir="outData/island/",gff="outData/prodigal/prodigal.gff",workdir=config["workdir"]
      
  shell: """
         rm -f {output.out}
         
         for i in `cat {input}|grep ">"|cut -c 2-`
         do
            cat {input}|seqkit grep -p ${{i}} > {params.outdir}/${{i}}.fasta
            cat {params.gff}|grep ${{i}} > {params.outdir}/${{i}}.gff
            python scripts/gff_to_gbk.py {params.outdir}/${{i}}.gff {params.outdir}/${{i}}.fasta
            mv {params.outdir}/${{i}}.gb {params.outdir}/${{i}}.gbk

            docker run \
              -it --rm \
              -v {params.workdir}/{params.outdir}/:/data  brinkmanlab/islandpath:1.0.0 /data/${{i}}.gbk /data/${{i}}_GIs.txt 

            cat {params.outdir}/${{i}}_GIs.txt | awk -v contig="$i" 'BEGIN{{FS=OFS="\t"}}{{print contig,$0}}' >> {output.out}

         done 

         cat {output.out}|awk 'BEGIN{{FS=OFS="\t"}}{{$2="Genomic Island";print $0}}' |sed "1 i Sequence\tFeature\tStart\tEnd" > {output.island}

         # gff提取gene坐标
         awk 'BEGIN{{FS=OFS="\t"}}
         NR==FNR&&FNR>1{{i++;sequence[i]=$1;start[i]=$3;end[i]=$4;next}}
         {{for (i in start) 
          {{if($4>=start[i]&&$5<=end[i]) print sequence[i],"Dimob","Gene",$4,$5,".",$7,$9,start[i],end[i],end[i]-start[i]+1}}}}' {params.outdir}/islands.txt {params.gff} | awk 'BEGIN {{OFS="\t"}}
          {{print $0}}'|sed "1 i #Sequence\tSource\tFeature\tGene_Start\tGene_End\tScore\tStrand\tGene_ID\tIsland_Start\tIsland_End\tIsland_Length" > {output.GeneIsland}


  """

rule Transposon:
  input: config["fasta"]
  output: "outData/transposon/Transposon.txt"
  params: outdir="outData/transposon/"
  shell: """
         mkdir -p {params.outdir}
         cd {params.outdir}

         /data1/software/TransposonPSI_08222010/transposonPSI.pl ../../{input} nuc

         sed "1 i Sequence\tSource\tMethod\tStart\tEnd\tScore\tStrand\tPhase\tAttributes" *TPSI.allHits.chains.gff3 > Transposon.txt
   

  """

rule casFinder:
  input: config["fasta"]
  output: "outData/crispr2/Crisprs.txt"
  params: outdir="outData/crispr/",workdir=config["workdir"]
  shell: """
           mkdir -p {params.outdir}
           docker run \
             -it --rm \
             -v {params.workdir}/:/data  icey_casfinder /bin/bash -c 'rm -rf /data/outData/crispr// && source ~/.profile && /opt/CRISPRCasFinder/CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def General -cas -i /data/{input} -out /data/{params.outdir} –keep -so /opt/CRISPRCasFinder/sel392v2.so'
           cp outData/crispr/TSV/Crisprs_REPORT.tsv outData/crispr2/Crisprs.txt

  """

rule prokka:
  input: config["fasta"]
  output: "outData/prokka/prokka.gbk"
  params:outdir="outData/prokka/"
  shell: """
           prokka {input} --outdir {params.outdir} --prefix prokka --force

  """  

rule prophage:
  input: gbk="outData/prokka/prokka.gbk"
  output: "outData/prophage/prophage.txt"
  params: outdir="outData/prophage/"
  conda:"../envs/R.yaml"
  shell: """

           /data1/software/miniconda2/bin/python /data1/software/PhiSpy/genbank_to_seed.py {input.gbk} {params.outdir}/organism_directory

           /data1/software/miniconda2/bin/python /data1/software/PhiSpy/PhiSpy.py -i {params.outdir}/organism_directory -o {params.outdir} 

           cat {params.outdir}/prophage_tbl.txt|awk 'BEGIN{{FS=OFS="\t"}}FNR==1{{print $0}}FNR>1{{if ($8==1)print $0}}' > {output}

  """  




