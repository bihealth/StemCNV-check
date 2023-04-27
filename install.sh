#!/usr/bin/bash

# TODO: this needs hook / conda set properly
# ## Conda env 
# mamba env create -f conda_env.yaml
# conda activate cnv-pipeline

extra_packages=$CONDA_PREFIX/pipeline
mkdir $extra_packages
start_path=$PWD

# this probably can not be added to bioconda
# Illumina Gencall software
cd $extra_packages
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
cd $CONDA_PREFIX/bin
ln -s $CONDA_PREFIX/pipeline/iaap-cli/iaap-cli
cd $start_path

#TODO
# bioconda now seems to have the bcftools 1.16 version available
# -> test that

# Build HTSlib & bcftools (incl gtc2vcf) from source. This is needed for gtc2vcf ( polysomy would be nice but doesn't work)
# Could also use bcftools from bioconda, but plugin is not up-do-date there & polysomy doesn't work either

## gtc2vcf // bcftools ##
# Build bcftools with gtc2vcf plugin from source
cd $extra_packages
git clone --branch=develop --recurse-submodules https://github.com/samtools/htslib.git
git clone --branch=develop https://github.com/samtools/bcftools.git
#Newest version broke this for me
#wget -P bcftools/plugins https://raw.githubusercontent.com/freeseek/gtc2vcf/master/{gtc2vcf.{c,h},affy2vcf.c}
# Alternative with last working commit fixed
git clone https://github.com/freeseek/gtc2vcf.git
cd gtc2vcf
git checkout 8c2f77469417dd8648d68b1249052fc9dd1bf44e
cd ..
cp gtc2vcf/gtc2vcf.c gtc2vcf/gtc2vcf.h gtc2vcf/affy2vcf.c bcftools/plugins/.

cd bcftools && autoheader && autoconf && ./configure --prefix=$CONDA_PREFIX && make && make install && cd ..
#cd bcftools && autoheader && autoconf && ./configure --enable-libgsl -prefix=$CONDA_PREFIX/bin && make && make install && cd ..
#cd bcftools && autoheader && autoconf && ./configure --enable-libgsl -prefix=$CONDA_PREFIX/bin && make clean && make USE_GPL=1 -fsanitize=address && make prefix=$CONDA_PREFIX/bin install && cd ..

cd bin
ln -s $CONDA_PREFIX/pipeline/bcftools/bcftools
cd $start_path

#TODO: maybe go ahead and add this to bioconda ?
## PennCNV - install algorithms from source
cd $extra_packages
wget https://github.com/WGLab/PennCNV/archive/v1.0.5.tar.gz
tar xvfz v1.0.5.tar.gz
#Makefile changes (add -Wno-error=format-security)
cp $start_path/PennCNV_overrides/Makefile PennCNV-1.0.5/kext/. 
cd PennCNV-1.0.5/kext && make
cd $CONDA_PREFIX/bin
ln -s $CONDA_PREFIX/pipeline/PennCNV-1.0.5/detect_cnv.pl PennCNV_detect
ln -s $CONDA_PREFIX/pipeline/PennCNV-1.0.5/kext
cd $start_path

# Another segmentation algorithm - specifically for array data
# R-package is too old / not found by conda, so use github
Rscript -e 'devtools::install_github("isglobal-brge/R-GADA")'
Rscript -e 'devtools::install_github("isglobal-brge/MAD")'

## Set Conda env vars
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d

echo "#!/bin/sh
export BCFTOOLS_PLUGINS=$extra_packages/bcftools/plugins" > ./etc/conda/activate.d/env_vars.sh

echo "#!/bin/sh
unset BCFTOOLS_PLUGINS" > ./etc/conda/deactivate.d/env_vars.sh

##close shop
cd $start_path
