# AUGUSTUS Training
# 将 GFF3 文件转换为 GeneBank 格式
ln -s ../../00.incipient_data/data_for_gene_prediction_and_RNA-seq/genome.fasta ./
gff2gbSmallDNA.pl ../../pasa/best_candidates.lowIdentity.gff3 genome.fasta 800 genes.raw.gb
# 去除错误的基因
new_species.pl --species=for_bad_genes_removing
etraining --species=for_bad_genes_removing --stopCodonExcludedFromCDS=false genes.raw.gb 2> train.err
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
filterGenes.pl badgenes.lst genes.raw.gb > genes.gb





# 将用于 training 的基因分成 2份。由于示例中的基因模型数量太少，因此使用来自GeneMark-ET对全基因组进行分析后得到的 GeneBank 文件。
rm genes.gb 
cp ~/00.incipient_data/data_for_gene_prediction_and_RNA-seq/genes.gb ./
randomSplit.pl genes.gb 100
new_species.pl --species=neurospora_crassa_OR74A
# 第一次training
etraining --species=neurospora_crassa_OR74A genes.gb.train > train.out
perl -e 'open IN, "train.out"; while (<IN>) { $tag = $1 if m/tag:.*\((.*)\)/; $taa = $1 if m/taa:.*\((.*)\)/; $tga = $1 if m/tga:.*\((.*)\)/; } while (<>) { s#/Constant/amberprob.*#/Constant/amberprob                   $tag#; s#/Constant/ochreprob.*#/Constant/ochreprob                   $taa#; s#/Constant/opalprob.*#/Constant/opalprob                    $tga#; print }' /opt/biosoft/augustus-3.2.2/config/species/neurospora_crassa_OR74A/neurospora_crassa_OR74A_parameters.cfg > 11
mv 11 /opt/biosoft/augustus-3.2.2/config/species/neurospora_crassa_OR74A/neurospora_crassa_OR74A_parameters.cfg
augustus --species=neurospora_crassa_OR74A genes.gb.test | tee firsttest.out
# 使用optimize_augustus.pl进行循环training找最优参数
# genes.gb.train中包含535个基因模型。对其再次进行分割，取其中 400 个基因模型用于对HMM模型参数进行优化时的准确性检测，剩下135个仅加入到training的过程中。本次优化过程如下：对某一个参数进行优化的时候，将400个基因模型随机分成8份，每份的基因模型数目为50个；取其中7份的基因模型和135个基因模型用于etraining，再使用剩下的1份基因模型进行准确性检测；总共有8份数据，这8份数据并行化运行，得到8个准确性检测值，取其均值用于参数的优化。
# 因此，为了能更准确快速的进行HMM参数文件优化，则需要注意两点：（1）根据计算资源确定并行数。当然，并行数越大，则准确性检测的值更准确。（2）用于准确性检测的基因模型数目。该数目除以并行数的值不要大于100，否则，每次准确性检测耗时会很长。当然，用于准确性检测的基因模型数目除以并行数的值越大，则准确性检测的值更准确。这个需要在准确性和耗时中权衡该值，推荐50-100。
# 若用于training的基因模型数目很多，比如有6000个；同时计算资源达160线程。则可以设置：将6000个基因模型分成1000和5000个；其中1000个用于参数文件的准确性检测；5000个用于training，设置并行化数100。
randomSplit.pl genes.gb.train 135
ln -s genes.gb.train.train training.gb.onlytrain
optimize_augustus.pl --species=neurospora_crassa_OR74A --rounds=5 --cpus=8 --kfold=8 --onlytrain=training.gb.onlytrain genes.gb.train.test > optimize.out
# 计算耗时~100min。
# accuracy value = (3*nucleotide_sensitivity + 2*nucleotide_specificity + 4*exon_sensitivity + 3*exon_specificity + 2*gene_sensitivity + 1*gene_specificity)/15
# Commonly observed values at this position range from 40 to 60 percent. If you obtain a very low value, this gives a strong indication that the obtained parameter set is not very useful for predicting genes accurately.
# 第二次training
etraining --species=neurospora_crassa_OR74A genes.gb.train
augustus --species=neurospora_crassa_OR74A genes.gb.test | tee secondtest.out
# CRF(Conditional Random Field) Training
cd /opt/biosoft/augustus-3.2.2/config/species/neurospora_crassa_OR74A
cp neurospora_crassa_OR74A_exon_probs.pbl neurospora_crassa_OR74A_exon_probs.pbl.withoutCRF
cp neurospora_crassa_OR74A_igenic_probs.pbl neurospora_crassa_OR74A_igenic_probs.pbl.withoutCRF
cp neurospora_crassa_OR74A_intron_probs.pbl neurospora_crassa_OR74A_intron_probs.pbl.withoutCRF
cd -
etraining --species=neurospora_crassa_OR74A --CRF=1 genes.gb.train
# 计算耗时~26min。
augustus --species=neurospora_crassa_OR74A genes.gb.test | tee secondtest.out.withCRF
# 比较 CRF 和 非 CRF 两种情况下的精确性。一般情况下，CRF training 的精确性要高些。若 CRF training 的精确性低些，则将备份的参数文件还原回去即可。
cd /opt/biosoft/augustus-3.2.2/config/species/neurospora_crassa_OR74A
cp neurospora_crassa_OR74A_exon_probs.pbl neurospora_crassa_OR74A_exon_probs.pbl.withCRF
cp neurospora_crassa_OR74A_igenic_probs.pbl neurospora_crassa_OR74A_igenic_probs.pbl.withCRF
cp neurospora_crassa_OR74A_intron_probs.pbl neurospora_crassa_OR74A_intron_probs.pbl.withCRF
cd -