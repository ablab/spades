#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run
# docker run -t -i -v <path with algorithmic-biology repository>:/gitrep --rm phusion/holy-build-box-64:latest  bash /gitrep/algorithmic-biology/assembler/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate

set -x

yum install bzip2-devel

cd /gitrep/algorithmic-biology/assembler
./spades_compile.sh -DSPADES_STATIC_BUILD=1

cd build_spades
make package

#ldd /gitrep/algorithmic-biology/assembler/bin/spades
#ldd /gitrep/algorithmic-biology/assembler/bin/hammer
#ldd /gitrep/algorithmic-biology/assembler/bin/bwa-spades
#ldd /gitrep/algorithmic-biology/assembler/bin/corrector
#ldd /gitrep/algorithmic-biology/assembler/bin/dipspades
#ldd /gitrep/algorithmic-biology/assembler/bin/ionhammer
#ldd /gitrep/algorithmic-biology/assembler/bin/scaffold_correction
