join_paired_ends.py -f ba9_S153_L001_R1_001.fastq.gz -r ba9_S153_L001_R2_001.fastq.gz  -o test_ba9_raw

cd test_ba9_raw

split_libraries_fastq.py  -i fastqjoin.join.fastq  -o split -q 19 --barcode_type 'not-barcoded'  --sample_ids ba9
