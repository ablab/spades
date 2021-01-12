#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run once:
# udocker load -i /home/akorobeynikov/spades-buildbox.tar
# udocker create --name=buildbox cab/spades-buildbox6-64:2.0.0-a4
# cd /tmp
# curl -LO https://github.com/Kitware/CMake/releases/download/v3.19.2/cmake-3.19.2-Linux-x86_64.tar.gz
# tar zxvf cmake-3.19.2-Linux-x86_64.tar.gz
# udocker run  --user=$(id -u) --volume=/tmp/:/temp buildbox cp -rv /temp/cmake-3.19.2-Linux-x86_64 /hbb/share/cmake-3.19
# Then run:
# udocker run  --user=$(id -u) --volume=<path to SPAdes repository>:/spades --cpuset-cpus=16 buildbox bash /spades/assembler/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate
export PATH=/hbb/share/cmake-3.19/bin/:$PATH

set -x

cd /spades/assembler
./spades_compile.sh -DSPADES_STATIC_BUILD=1

cd build_spades
make package
