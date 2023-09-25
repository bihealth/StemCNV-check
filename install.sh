#!/usr/bin/bash

# TODO: this needs hook / conda set properly
# ## Conda env 
# mamba env create -f conda_env.yaml
# conda activate cnv-pipeline

extra_packages=$CONDA_PREFIX/pipeline
mkdir $extra_packages
start_path=$PWD
# If an http proxy is set we need it to apply to ftp as well
if [[ -n $http_proxy && -z ${ftp_proxy+x} ]]; then
  export ftp_proxy=$http_proxy
fi


# Illumina Gencall software
# this probably can not be added to bioconda
cd $extra_packages
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
cd $CONDA_PREFIX/bin
ln -s ../pipeline/iaap-cli/iaap-cli
cd $start_path

## PennCNV - install algorithms from source
#TODO: maybe go ahead and add PennCNV to bioconda ?
# fix needed since this was build for gcc <= 2.26
# -> can probably be changed by adapting khmm (.h / .c / .o) or khmm_wrap file

# perl issues
# 5.32 -> gcc can't find xlocale whn including perl c code into PennCNV
# but works if that is mitigated:
sed -i 's/xlocale.h/locale.h/' $CONDA_PREFIX/lib/perl5/5.32/core_perl/CORE/perl.h
# all other perl versions -> some sort of ld // perl version problem when building PennCNV
#/usr/bin/ld: [...]/libperl.a(op.o): relocation R_X86_64_32 against `.rodata.str1.8' can not be used when making a shared object; recompile with -fPIC
#/usr/bin/ld: failed to set dynamic section sizes: bad value
#collect2: error: ld returned 1 exit status

cd $extra_packages
wget https://github.com/WGLab/PennCNV/archive/v1.0.5.tar.gz
tar xvfz v1.0.5.tar.gz
#Makefile changes (add -Wno-error=format-security)
cp $start_path/PennCNV_overrides/Makefile PennCNV-1.0.5/kext/. 
cd PennCNV-1.0.5/kext && make
cd $CONDA_PREFIX/bin
ln -s ../pipeline/PennCNV-1.0.5/detect_cnv.pl PennCNV_detect
ln -s ../pipeline/PennCNV-1.0.5/kext
cd $start_path

# Likely needed for pdf generation on WSL - doesn't take long though so better to include
Rscript -e 'tinytex::install_tinytex()'

#TODO - keep these?

## R-GADA
# Another segmentation algorithm - specifically for array data
# R-package in conda is too old / not found by conda, so use github
Rscript -e 'devtools::install_github("isglobal-brge/R-GADA")'
Rscript -e 'devtools::install_github("isglobal-brge/MAD")'

# Not yet available in conda
pip install 'gosling[all]'

##close shop
cd $start_path
