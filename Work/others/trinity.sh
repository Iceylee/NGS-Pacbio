########################### Trinity Comstruction #############################
##############################################################################
### Trinity 
/opt/trinityrnaseq-Trinity-v2.4.0/Trinity --seqType fq --max_memory 250G --CPU 30 --samples_file sample_list.txt  --output /data1/ClientData/2018_03/HuGaoSheng/1.Trinity

### N50 
/opt/trinityrnaseq-Trinity-v2.4.0/util/TrinityStats.pl Trinity.fasta

### CD-HIT 0.9
mkdir Trinity_CD-HIT_0.9;cd Trinity_CD-HIT_0.9
cd-hit-est -i ../Trinity.fasta -o Trinity_CD-HIT_0.9.fa  -n 8 -p 1 -g 1 -M 0 -T 20 -d 0

########################### Trinity Mapping ##################################
##############################################################################
### prepase Index  Trinity-CD-HIT.fa 
/opt/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --transcripts Trinity_CD-HIT_0.9.fa  --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

for i in $(seq 180315 180323)
do
        /opt/trinityrnaseq-Trinity-v2.4.0/util/align_and_estimate_abundance.pl --thread_count 20 --transcripts /data1/ClientData/2018_03/HuGaoSheng/1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa  --seqType fq --left "/DataBackUp/2018_03/HuGaoSheng/clean_data/WR"$i"S_clean_R1.fastq.gz" --right "/DataBackUp/2018_03/HuGaoSheng/clean_data/WR"$i"S_clean_R2.fastq.gz" --est_method RSEM --aln_method bowtie --trinity_mode --SS_lib_type RF --gene_trans_map WR$iS --output_dir  "/data1/ClientData/2018_03/HuGaoSheng/2.RSEM/WR"$i"S_rsem_outdir"

done


############ Trinity Find Coding Regions Within Transcripts ##################
##############################################################################
### TransDecoder
/data/software/TransDecoder/TransDecoder.LongOrfs  -t ../Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa

hmmscan --cpu 8 --domtblout pfam.domtblout /data/Pfam/Pfam-A.hmm longest_orfs.pep

/data/software/TransDecoder/TransDecoder.Predict -t ../Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --retain_pfam_hits  Trinity_CD-HIT_0.9.fa.transdecoder_dir/pfam.domtblout --cpu 20


###################### Trinity Gene Annotation ###############################
##############################################################################
## Uniprot
/data/software/diamond-linux64/diamond makedb --in uniprot_sprot.fasta -d uniprot_sprot

/data/software/diamond-linux64/diamond blastx -q ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa -d /data/Uniprot/uniprot_sprot -o Uniprot_sport.out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

## KOG
/data/software/diamond-linux64/diamond blastx -q ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa -d /data/KOGDatabase/clean.kog  -o KOG_annotation.bk.out  -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle

### Diamond blastx ### 
/data/software/diamond-linux64/diamond blastx -q Trinity_CD-HIT_0.9.fa -d  /data/DiamondDB/nr -o ../../3.Annotation/NR/Trinity_CD-HIT_0.9_blastx_nr.txt -p 20 -e 1e-5 --max-target-seqs 1 --more-sensitive --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle




### emapper ###  
/data/software/eggnog-mapper/emapper.py -i ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --output Trinity_CD-HIT_0.9_maNOG -m diamond --cpu 20


####
#  emapper-1.0.3
# ./emapper.py  -i ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --output Trinity_CD-HIT_0.9_maNOG -m diamond --cpu 20
  /data/software/eggnog-mapper/bin/diamond blastp -d /data/software/eggnog-mapper/data/eggnog_proteins.dmnd -q /data1/ClientData/2018_03/HuGaoSheng/1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --more-sensitive --threads 20 -e 0.001000 -o /data1/ClientData/2018_03/HuGaoSheng/3.Annotation/emappertmp_dmdn_Ysb_RS/c7d80423f3a44f0a930fb44a51f4dacd --top 3
####

### eggnog blastp 
/data/software/eggnog-mapper/bin/diamond blastp -d /data/software/eggnog-mapper/data/eggnog_proteins.dmnd -q /data1/ClientData/2018_03/HuGaoSheng/1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --more-sensitive --threads 20 -e 0.001000 -o /data1/ClientData/2018_03/HuGaoSheng/3.Annotation/emapper/blastp_Trinity_CD-hit.out --top 3 
####




###################### Trinity Gene Expression ###############################
##############################################################################
### Gene
/opt/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix RSEM --name_sample_by_basedir WR180315S_rsem_outdir/RSEM.genes.results WR180316S_rsem_outdir/RSEM.genes.results WR180317S_rsem_outdir/RSEM.genes.results WR180318S_rsem_outdir/RSEM.genes.results WR180319S_rsem_outdir/RSEM.genes.results WR180320S_rsem_outdir/RSEM.genes.results WR180321S_rsem_outdir/RSEM.genes.results WR180322S_rsem_outdir/RSEM.genes.results WR180323S_rsem_outdir/RSEM.genes.results

### Isoform
/opt/trinityrnaseq-Trinity-v2.4.0/util/abundance_estimates_to_matrix.pl --est_method RSEM  --out_prefix RSEM_isoform --name_sample_by_basedir ../WR180315S_rsem_outdir/RSEM.isoforms.results ../WR180316S_rsem_outdir/RSEM.isoforms.results ../WR180317S_rsem_outdir/RSEM.isoforms.results ../WR180318S_rsem_outdir/RSEM.isoforms.results ../WR180319S_rsem_outdir/RSEM.isoforms.results ../WR180320S_rsem_outdir/RSEM.isoforms.results ../WR180321S_rsem_outdir/RSEM.isoforms.results ../WR180322S_rsem_outdir/RSEM.isoforms.results ../WR180323S_rsem_outdir/RSEM.isoforms.results



###################### Trinity Diff Expression ###############################
##############################################################################
###DESeq2 ### Diff Analysis use Gene ### If you want see isoform,just use "RSEM_isoform.counts.matrix" 
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix RSEM.counts.matrix --method DESeq2 --samples_file sample_list.txt --output DESeq2_results

/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../RSEM.TMM.EXPR.matrix --samples ../sample_list.txt -P 1e-3 -C 2 --output cluster_results 


################
TMM to FPKM


cut -f 1,3,4 sampleA.RSEM.isoforms.results > feature_lengths.txt
$TRINITY_HOME/Analysis/DifferentialExpression/run_TMM_normalization_write_FPKM_matrix.pl --matrix counts.matrix --lengths feature_lengths.txt

################




### Cluster ###
 /opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/define_clusters_by_cutting_tree.pl -K 6 -R cluster_results.matrix.RData 




### Filter isoform TPM ###
/opt/trinityrnaseq-Trinity-v2.4.0/util/filter_low_expr_transcripts.pl --matrix ../2.RSEM/Matrix_isoforms/RSEM_isoform.TPM.not_cross_norm  --transcripts ../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa --min_expr_any 10.0    > Trinity.filtered.fa 
#### Retained 30763 / 618201 = 4.98% of total transcripts. ####



### Running RNAMMER to identify rRNA transcripts
/data/software/Trinotate-Trinotate-v3.1.1/util/rnammer_support/RnammerTranscriptome.pl --transcriptome Trinity.filtered.fa --path_to_rnammer /data/software/RNAMMER/rnammer


###################### Diff Expression  KEGG  ################################
##############################################################################
#### 
python GetDiffGeneAnnotation.py RSEM.counts.matrix.Cambium_vs_ScaleLeaf.DESeq2.DE_results.P1e-3_C2.DE.subset ../../../3.Annotation/KEGG/Trinity.filtered_KAAS_annotation_gene.txt > Cambium_vs_ScaleLeaf.DiffGene.KEGG.txt 


### Get Seq Length
/opt/trinityrnaseq-Trinity-v2.4.0/util/misc/fasta_seq_length.pl ../../1.Trinity/Trinity_CD-HIT_0.9/Trinity_CD-HIT_0.9.fa  > Trinity_CD-HIT_0.9.seqLengths

### PCA
/opt/trinityrnaseq-Trinity-v2.4.0/Analysis/DifferentialExpression/PtR --matrix RSEM.counts.matrix   -s ../sample_list.txt  --sample_cor_matrix  --prin_comp 3
