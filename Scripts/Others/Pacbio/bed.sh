git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make clean; make all

# create a BED file with the region of interest
# the spaces in the BED file are tabs
cat my_region.bed

# calculate the coverage in this region
bedtools coverage -a my_region.bed -b wgEncodeRikenCageHchCellPapAlnRep1.bam -bed

# calculate the coverage per base pair in this region
# on the same strand
bedtools coverage -a my_region.bed -b wgEncodeRikenCageHchCellPapAlnRep1.bam -bed -d -s | gzip > wgEncodeRikenCageHchCellPapAlnRep1_chr22_my_region_pos.tsv.gz
# on the opposite strand
bedtools coverage -a my_region.bed -b wgEncodeRikenCageHchCellPapAlnRep1.bam -bed -d -S | gzip > wgEncodeRikenCageHchCellPapAlnRep1_chr22_my_region_neg.tsv.gz
 
# check out the files
gunzip -c wgEncodeRikenCageHchCellPapAlnRep1_chr22_my_region_pos.tsv.gz | head -5

gunzip -c wgEncodeRikenCageHchCellPapAlnRep1_chr22_my_region_neg.tsv.gz | head -5