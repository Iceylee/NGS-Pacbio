
rule circos_data:
  input:rRNA="outData/rRNA/rRNA.gff",tRNA="outData/tRNA/tRNA.gff",sRNA="outData/sRNA/sRNA.gff",gff="outData/prodigal/genome.gff",fasta=config["fasta"],crispr="outData/crispr2/Crisprs.txt",island="outData/island/islands.txt",prophage="outData/prophage/prophage.txt",color=config["color"]
  output:"outData/circos/data/"
  params:datadir="outData/circos/data/"
  shell:"""
        seqkit fx2tab -l -g -n -i -H {input.fasta} |cut -f 1,4|sed '1 d'|sort -V -k1,1 > {params.datadir}/temp1
        awk 'FNR==NR{{i=i+1;a[i]=$1;next}}{{i=FNR;split($1,x,"_");print "chr - "$1" "x[2]" 0 "$2" "a[i]}}' {input.color} {params.datadir}/temp1> {params.datadir}/karyotype.txt
        
        GC_percent=`cat {input.gff}|sed -n 3p|grep -E "gc_cont=[0-9]+\.[0-9]+" -o|grep -E -o "[0-9]+\.[0-9]+"`
        GC=`echo "$GC_percent * 0.0100"|bc`

        python scripts/GCcalc.py -f {input.fasta} -w 10000 -s 10000 >{params.datadir}/gc_content_skew.txt
        cat {params.datadir}/gc_content_skew.txt|awk -v var="$GC" '{{FS=OFS="\t"}}{{print $1,$2,$3,$4-var}}'>{params.datadir}/gc_count.txt
        cat {params.datadir}/gc_content_skew.txt|awk '{{FS=OFS="\t"}}{{print $1,$2,$3,$5}}'>{params.datadir}/gc_skew.txt

        cat {input.gff} |awk '$7=="+"&&FNR>1{{FS=OFS="\t";print $1,$4,$5}}' >{params.datadir}/sensegene.gff
        cat {input.gff} |awk '$7=="-"&&FNR>1{{FS=OFS="\t";print $1,$4,$5}}' >{params.datadir}/antigene.gff


        awk 'NR==FNR{{i=i+1;a[i]=$1;next}}{{num=int(rand()*7)+1;print $1,$2,$3,"fill_color="a[num]"_a2"}}' scripts/color2.list {params.datadir}/sensegene.gff >{params.datadir}/sensegene_color.gff
        awk 'NR==FNR{{i=i+1;a[i]=$1;next}}{{num=int(rand()*7)+8;print $1,$2,$3,"fill_color="a[num]"_a2"}}' scripts/color2.list {params.datadir}/antigene.gff >{params.datadir}/antigene_color.gff

        cat {input.rRNA}|grep -v "Sequence"|cut -f 1,4,5 >{params.datadir}/rRNA_plot.gff
        cat {input.tRNA}|grep -v "Sequence"|cut -f 1,4,5 >{params.datadir}/tRNA_plot.gff
        cat {input.sRNA}|grep -v "Sequence"|cut -f 1,4,5 >{params.datadir}/sRNA_plot.gff

        cat {input.prophage}|sed '1d'|cut -f 3,4,5 > {params.datadir}/Prophage_plot.txt
        cat {input.crispr}|sed '1d'|cut -f 2,6,7 > {params.datadir}/CRISPR_plot.txt
        cat {input.island}|sed '1d'|cut -f 1,3,4 > {params.datadir}/GIs_plot.txt

  """

rule circos:
  input:"outData/circos/data/"
  output:datadir="outData/circos/",plot="outData/circos/circos.png"
  params:contig=config["contig"],short_chr=config["short_chr"]
  shell:"""
  		if [ {params.contig} = "TRUE" ]; then
  			cp -r scripts/one_contig/plot {output.datadir}/
  		else
  			cp -r scripts/multiple_contigs/plot {output.datadir}/

  			#short chr draw no ticks
  			sed -i '7ichromosomes={params.short_chr}' {output.datadir}/plot/ticks.conf
  		fi

  		cd {output.datadir}
        
  		circos -conf plot/circos.conf 
  	"""









