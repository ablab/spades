#!/bin/bash

# Script for compiling SPAdes using holy build box
# To compile spades run once:
# udocker load -i /Sid/prjbel/ablab.holy-build-box-x64.latest.tar
# udocker create --name=ablab-buildbox ghcr.io/ablab/holy-build-box-x64
# Then run:
# udocker run  --user=$(id -u) --volume=<path to SPAdes repository>:/spades --cpuset-cpus=16 buildbox bash /spades/assembler/src/tools/spades_holybox_compile.sh

set -e
source /hbb_exe/activate
export PATH=/hbb/share/cmake-3.19/bin/:$PATH

set -x

cd /spades
#./spades_compile.sh -DSPADES_STATIC_BUILD=1
mkdir -p build_spades
cd build_spades
cmake ../src -DSPADES_STATIC_BUILD=ON -DZLIB_ROOT=/hbb_exe -DOpenMP_pthread_LIBRARY=“-lpthread” -DBZIP2_LIBRARIES=/hbb_exe/lib/libbz2.a -DBZIP2_INCLUDE_DIR=/hbb_exe/include
make -j 8
make package
