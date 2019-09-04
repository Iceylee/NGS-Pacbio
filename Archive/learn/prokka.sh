prokka assembly.fasta

#prodigal  (CDS)
#aragorn(tRNA) Barrnap/RNAmmer(rRNA) 
#ncRNA(rfam) 需参数 --rfam
#CRISPR 
#anno (uniprot??)


#宏基因组
prokka final.contigs.fa --outdir prokka_annotation --prefix metagG --metagenome --kingdom Bacteria

#
prokka 

python /data1/software/PhiSpy/genbank_to_seed.py PROKKA_09142018.gbk organism_directory
/data1/software/PhiSpy/PhiSpy.py -i organism_directory -o output_directory -c









Beginner

# Vanilla (but with free toppings)
% prokka contigs.fa

# Look for a folder called PROKKA_yyyymmdd (today's date) and look at stats
% cat PROKKA_yyyymmdd/*.txt
Moderate

# Choose the names of the output files
% prokka --outdir mydir --prefix mygenome contigs.fa

# Visualize it in Artemis
% art mydir/mygenome.gff
Specialist

# Have curated genomes I want to use to annotate from
% prokka --proteins MG1655.gbk --outdir mutant --prefix K12_mut contigs.fa

# Look at tabular features
% less -S mutant/K12_mut.tsv
Expert

# It's not just for bacteria, people
% prokka --kingdom Archaea --outdir mydir --genus Pyrococcus --locustag PYCC

# Search for your favourite gene
% exonerate --bestn 1 zetatoxin.fasta mydir/PYCC_06072012.faa | less
Wizard

# Watch and learn
% prokka --outdir mydir --locustag EHEC --proteins NewToxins.faa --evalue 0.001 --gram neg --addgenes contigs.fa

# Check to see if anything went really wrong
% less mydir/EHEC_06072012.err

# Add final details using Sequin
% sequin mydir/EHEC_0607201.sqn