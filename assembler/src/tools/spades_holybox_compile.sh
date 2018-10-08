#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run
# docker run -u $(id -u):$(id -g) -t -i -v <path to SPAdes repository:/spades/ --rm cab/spades-buildbox6-64:2.0.0-a4  bash /spades/assembler/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate

set -x

cd /spades/assembler
./spades_compile.sh -DSPADES_STATIC_BUILD=1

cd build_spades
make package
