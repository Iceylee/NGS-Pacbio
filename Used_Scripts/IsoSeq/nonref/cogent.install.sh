conda create -n anaCogent python=2.7 anaconda
source activate anaCogent
conda install -n anaCogent biopython
conda install -n anaCogent -c http://conda.anaconda.org/cgat bx-python
#conda install -c http://conda.anaconda.org/timcera pulp
conda install -c agerlach pulp

#Icey
pip install networkx
pip install scipy
pip install scikit-image #matplotlib 一起装了


#
cd <your_dir>
git clone https://github.com/Magdoll/Cogent.git
cd Cogent
git submodule update --init --recursive
cd  Complete-Striped-Smith-Waterman-Library/src
make
export LD_LIBRARY_PATH=~/biosoft/Cogent/Complete-Striped-Smith-Waterman-Library/src
export PYTHONPATH=~/biosoft/Cogent/Complete-Striped-Smith-Waterman-Library/src
cd ../../
python setup.py build
python setup.py install #warning
#test
run_mash.py --version
run_mash.py Cogent 1.0
python setup.py test #test_cogent.py报错 
python ~/biosoft/Cogent/Cogent/test/test_cogent.py #ok


#报错error: pyparsing 1.5.7 is installed but pyparsing!=2.0.4,!=2.1.2,!=2.1.6,>=2.0.1 is required by set(['matplotlib'])
pkg_resources.DistributionNotFound: The 'pyparsing<=1.9.9' distribution was not found and is required by pulp

pip install pyparsing==1.5.7


##diamond
wget http://github.com/bbuchfink/diamond/releases/download/v0.9.18/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
mv diamond bin/
#索引
nr uniprot_trembl.fasta Pfam-A-duplicate.fasta(路径)
nohup diamond makedb --in Pfam-A-duplicate.fasta -d Pfam-A-duplicate -p 22 &


##
diamond blastx -d nr -q reads.fna -o nr.m8 --max-target-seqs 1 --evalue 1e-5 -p 20



