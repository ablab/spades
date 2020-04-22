#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run
# udocker load -i /home/akorobeynikov/spades-buildbox.tar
# udocker create --name=buildbox cab/spades-buildbox6-64:2.0.0-a4
# udocker run  --user=$(id -u) --volume=<path to SPAdes repository>:/spades --cpuset-cpus=16 buildbox bash /spades/assembler/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate

set -x

cd /spades/assembler
./spades_compile.sh -DSPADES_STATIC_BUILD=1

cd build_spades
make package
