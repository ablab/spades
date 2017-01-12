#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run
# docker run -u $(id -u):$(id -g) -t -i -v <path to spades_compile.sh>:/spades/ --rm cab/spades-buildbox  bash /spades/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate

set -x

cd /spades
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
