#!/usr/bin/bash

# This assumes to run within the scnv-quac conda environment
if [[ ! $CONDA_PREFIX == *scnv-quac ]]; then
  echo "Please run this script within the scnv-quac conda environment"
  exit 1
fi

extra_packages=$CONDA_PREFIX/pipeline
mkdir $extra_packages
start_path=$PWD
script_path=$(dirname -- "$( readlink -f -- "$0"; )";)
# If an http proxy is set we need it to apply to ftp as well
if [[ -n $http_proxy && -z ${ftp_proxy+x} ]]; then
  export ftp_proxy=$http_proxy
fi

## get Illumina Gencall software
cd $extra_packages
wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/iaap/iaap-cli-linux-x64-1.1.0.tar.gz
tar xzvf iaap-cli-linux-x64-1.1.0.tar.gz iaap-cli-linux-x64-1.1.0/iaap-cli --strip-components=1
cd $CONDA_PREFIX/bin
ln -s ../pipeline/iaap-cli/iaap-cli
cd $start_path

## PennCNV - build from source

# perl issues:
# - build with any perl version except 5.32 fail on some sort of "ld error"
# - building with perl 5.32 still needs this fix (sicne gcc can't find xlocale when including perl c code into PennCNV)
sed -i 's/xlocale.h/locale.h/' $CONDA_PREFIX/lib/perl5/5.32/core_perl/CORE/perl.h

cd $extra_packages
wget https://github.com/WGLab/PennCNV/archive/v1.0.5.tar.gz
tar xvfz v1.0.5.tar.gz
#Makefile changes (add -Wno-error=format-security)
cp $script_path/supplemental-files/penncnv-makefile PennCNV-1.0.5/kext/Makefile
cd PennCNV-1.0.5/kext && make
cd $CONDA_PREFIX/bin
ln -s ../pipeline/PennCNV-1.0.5/detect_cnv.pl PennCNV_detect
ln -s ../pipeline/PennCNV-1.0.5/kext
cd $start_path

##close shop
cd $start_path
