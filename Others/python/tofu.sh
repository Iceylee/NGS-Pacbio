
conda create -n anapython2.7.9 python=2.7.9
source activate anapython2.7.9

conda install virtualenv
export VENV_TOFU=~/biosoft/VENV_TOFU
virtualenv --system-site-packages -p python $VENV_TOFU
source $VENV_TOFU/bin/activate


conda install biopython
pip install --upgrade setuptools
pip install Cython
pip install numpy
pip install bx-python ###
pip install pysam
pip install pbcore

#fusionfinder
pip install networkx
pip install scipy

cd cDNA_primer-master/pbtranscript-tofu/pbtranscript/pbtools/pbtranscript/branch/C/
python setup2.py build_ext --inplace
cd ../../../../
make

touch /home/liyb/biosoft/VENV_TOFU/lib/python2.7/site-packages/pbtools.pbtranscript-2.2.3-py2.7-linux-x86_64.egg/pbtools/pbtranscript/modified_bx_intervals/__init__.py

cd ~/biosoft/cDNA_primer-master/pbtranscript-tofu/external_daligner
unzip DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod.zip
cd DALIGNER-d4aa4871122b35ac92e2cc13d9b1b1e9c5b5dc5c-ICEmod/
make
cp HPCdaligner HPCmapper LA4Ice LAcat LAcheck LAmerge LAshow LAsort LAsplit daligner ${VENV_TOFU}/bin

cd ../
unzip DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa.zip
cd DAZZ_DB-40bb7e4b2041cdfd3c50b22201301b8df06342fa/
make  #有几个报错 注意
cp Catrack DAM2fasta DB2fasta DB2quiva DBdust DBrm DBshow DBsplit DBstats fasta2DAM fasta2DB quiva2DB simulator ${VENV_TOFU}/bin