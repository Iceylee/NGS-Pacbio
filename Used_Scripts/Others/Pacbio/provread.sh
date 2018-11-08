# small genome (10Mbp) or small memory machine or short runtime limit
SeqChunker -s 20M -o pb-%03d.fq pb-subreads.fq

# large genome (>500Mbp) and >8 GB RAM, no runtime limit
# either go for one instance of proovread per SMRT cell or
SeqChunker -s 1G -o pb-%03d.fq pb-subreads.fq

proovread -l pb-001.fq -s reads.fq [-u unitigs.fa] --pre pb-001