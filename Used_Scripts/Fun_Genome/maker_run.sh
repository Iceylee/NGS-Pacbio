#获取gene-mark-es的hmm文件
#es:DNA序列直接预测
#et：结合转录组bed等文件来预测。更复杂
mkdir genemark_es_et
cd genemark_es_et
gmes_petap.pl --sequence ../genome.fasta --ES --fungus --cores 15
~/biosoft/PASApipeline-master/misc_utilities/gtf_to_gff3_format.pl genemark.gtf ../genome.fasta > genemark.gff3

#maker
# 准备输入文件
ln -s ../../00.incipient_data/data_for_gene_prediction_and_RNA-seq/genome.fasta ./
ln -s ../../00.incipient_data/data_for_gene_prediction_and_RNA-seq/Trinity-GG.fasta ./
ln -s ../../00.incipient_data/data_for_gene_prediction_and_RNA-seq/Cryphonectria_parasitica.fasta proteins.fasta
ln -s ../../00.incipient_data/data_for_gene_prediction_and_RNA-seq/RepeatModeler_database.fa.classified ./
ln -s /home/liyb/analysis/learn/maker/Chen/genemark_es_et/output/gmhmm.mod ./
ln -s ../../04.genome_feature_analysis/RNAmmer/Ncrassa/rRNA.fasta ./
# 准备配置文件
maker -CTL
perl -p -i -e 's#RepeatMasker=.*#RepeatMasker=/opt/biosoft/RepeatMasker/RepeatMasker#' maker_exe.ctl

#maker_opts.ctl
genome=genome.fasta
est=Trinity-GG.fasta
protein=proteins.fasta
model_org=fungi
rmlib=RepeatModeler_database.fa.classified
gmhmm=gmhmm.mod
augustus_species=neurospora_crassa_OR74A
est2genome=1
trna=1
correct_est_fusion=1
keep_preds=1

#
perl -p -i -e 's/^augustus_species=.*/augustus_species=chaetomium_globosum/' maker_opts.ctl

# 运行maker - 15cpu 3:55
maker

# 计算耗时~62min。
cd genome.maker.output
gff3_merge -d genome_master_datastore_index.log 
grep -P "\tmaker\t" genome.all.gff > genome.maker.gff3
maker_gff3_out_sort_rename.pl --maxAED 0.1 genome.maker.gff3 > genome.maker.AED_0.1.gff3