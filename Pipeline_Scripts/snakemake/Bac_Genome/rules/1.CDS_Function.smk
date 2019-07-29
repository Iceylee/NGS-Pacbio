
rule prodigal:
  input: config["fasta"]
  output: pep="outData/prodigal/prodigal.pep",
          cds="outData/prodigal/prodigal.cds",
          gff="outData/prodigal/prodigal.gff",
          stat="outData/prodigal/prodigal.stat"       
  message: """--- Prodigal."""
  shell: """
         prodigal -i {input} -d {output.cds} -f gff -o {output.gff} -s {output.stat} -a {output.pep}
  """

rule rename:
  input:pep="outData/prodigal/prodigal.pep",
        cds="outData/prodigal/prodigal.cds",
        gff="outData/prodigal/prodigal.gff"
  output:pep="outData/prodigal/genome.pep",
         cds="outData/prodigal/genome.cds",
         gff="outData/prodigal/genome.gff"
  params: name=config["name"]
  shell:"""

          awk -v name={params.name} -f scripts/pep2.awk {input.pep} >{output.pep}
          awk -v name={params.name} -f scripts/pep2.awk {input.cds} >{output.cds}
          awk -v name={params.name} -f scripts/gff2.awk {input.gff} >{output.gff} 
      
  """

DBS = config["db"]

rule diamond:
  input:pep="outData/prodigal/genome.pep"  
  output:
         expand("outData/diamond/{db}.anno.out",
               db=DBS)
  run:
        for i in DBS:
        
          cmd = "/data1/software/diamond-linux64/diamond blastp -q %s -d %s -o outData/diamond/%s.anno.out -p %s -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle" % (input, config["db"][i], i, config["threads"])

          os.system(cmd)
        

rule cog_count:
  input:"outData/diamond/cog.anno.out"
  output:"outData/diamond/cog_count_fun.txt"
  params:csv=config["cog"]["csv"],cognames=config["cog"]["cognames"],
         fun=config["cog"]["fun"],outdir="outData/diamond/"
  shell:"""
        bash scripts/cog.sh {params.csv} {params.cognames} {params.fun} {input} {output} {params.outdir}

  """

rule cog_plot:
  input:"outData/diamond/cog_count_fun.txt"
  output:barplot="outData/plots/COG_barplot.pdf"
  params:prefix="outData/plots/COG_barplot"
  conda:"../envs/R.yaml"
  shell:"""
        Rscript scripts/cog_plot.R {input} {params.prefix}

  """


rule go_count:
  input:"outData/diamond/sp.anno.out"
  output:anno="outData/diamond/go.anno.out",count="outData/diamond/go_count.out"
  params:tab=config["go"]["tab"],gaf=config["go"]["gaf"],outdir="outData/diamond/"
  shell:"""

        bash scripts/sp2go.sh {input} {params.tab} {params.gaf} {output.anno} {output.count} {params.outdir}       

  """

rule go_plot:
  input:"outData/diamond/go_count.out"
  output:barplot="outData/plots/GO_barplot.pdf"
  params:prefix="outData/plots/GO_barplot"
  conda:"../envs/R.yaml"
  shell:"""
        Rscript scripts/GO_plot.R {input} {params.prefix}

  """

rule kegg_count:
  input:config["kegg"]
  output:anno="outData/diamond/kegg.anno.out",count="outData/diamond/kegg_count.out"
  shell:"""

        bash scripts/kegg.sh {input} {output.anno} {output.count}     

  """


rule kegg_plot:
  input:"outData/diamond/kegg_count.out"
  output:barplot="outData/plots/KEGG_barplot.pdf"
  params:prefix="outData/plots/KEGG_barplot"
  conda:"../envs/R.yaml"
  shell:"""
        Rscript scripts/KEGG_plot.R {input} {params.prefix}

  """
















rule anno_summary:
  input:"outData/prodigal/genome.pep"
  output:summary="outData/anno/Annotation_Summary.txt",list="outData/anno/genome.list"
  params:indir="outData/diamond/",outdir="outData/anno/"
  conda:"../envs/R.yaml"
  shell:"""

        cat {input} |grep ">"|sed 's/>//g' |awk 'BEGIN{{FS=" "}}{{print $1}}' > {output.list}

        Rscript scripts/anno_summary.R {params.indir} {output.list} {params.outdir} 

  """
