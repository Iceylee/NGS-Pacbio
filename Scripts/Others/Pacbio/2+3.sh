#install  isoaux
# clone repository and its submodules
git clone --recursive https://github.com/bowhan/isoaux.git

# install
cd isoaux && mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE .. && make && make install
# if cmake is not available, install brew (see below) and install cmake with brew install cmake

# for rest of the third-party software, we recommend to use brew
# to install brew
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/install/master/install)"
# then update ${PATH}

# install 3rd party software with brew and pip
brew tap homebrew/science
brew install bwa bowtie2 gmap-gsnap seqtk bedtools samtools # rna-star
# Note: the current version of STAR_2.4.1d on science/brew is outdated
# Please obtain the current version (2.5) from source: https://github.com/alexdobin/STAR/releases

pip install RSeQC