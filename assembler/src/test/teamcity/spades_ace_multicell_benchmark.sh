#!/bin/bash

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# Copyright (c) 2011-2014 Saint Petersburg Academic University
# All Rights Reserved
# See file LICENSE for details.
############################################################################


./prepare_cfg
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "prepare_cfg finished with exit code $errlvl"
    exit $errlvl
fi

./spades_compile.sh
errlvl=$?
if [ "$errlvl" -ne 0 ]
then
    echo "spades_compile finished with exit code $errlvl"
    exit $errlvl
fi

timestamp="`date +%Y-%m-%d_%H-%M-%S`"
if [ "$#" > 0 ]
then
    timestamp=$timestamp"_"$1
fi

mkdir -p ./ace_benchmark

echo "Starting E.coli is220 on ace..."
rm -f ./ace_benchmark/ECOLI_IS220_QUAKE.log
./src/test/teamcity/teamcity.py /Johnny/data/input/Bacteria/E.coli/K12/is220/ECOLI_IS220_QUAKE_BM.info $1 >> ./ace_benchmark/ECOLI_IS220_QUAKE_$timestamp.log 2>> ./ace_benchmark/ECOLI_IS220_QUAKE_$timestamp.log
echo "finished with code "$?

echo "Starting E.coli is480 on ace..."
rm -f ./ace_benchmark/ECOLI_IS480_QUAKE.log
./src/test/teamcity/teamcity.py /Johnny/data/input/Bacteria/E.coli/K12/is480/ECOLI_IS480_QUAKE_BM.info >> ./ace_benchmark/ECOLI_IS480_QUAKE_$timestamp.log 2>> ./ace_benchmark/ECOLI_IS480_QUAKE_$timestamp.log 
echo "finished with code "$?

echo "Starting B.faecium on ace..."
rm -f ./ace_benchmark/BFAECIUM_QAUKE.log
./src/test/teamcity/teamcity.py /Johnny/data/input/Bacteria/B.faecium/BFAECIUM_QAUKE_BM.info >> ./ace_benchmark/BFAECIUM_QAUKE_$timestamp.log 2>> ./ace_benchmark/BFAECIUM_QAUKE_$timestamp.log 
echo "finished with code "$?

echo "Starting L.gasseri on ace..."
rm -f ./ace_benchmark/LGASSERI_QUAKE.log
./src/test/teamcity/teamcity.py /Johnny/data/input/Bacteria/L.gasseri/LGASSERI_QUAKE_BM.info >> ./ace_benchmark/LGASSERI_QUAKE_$timestamp.log 2>> ./ace_benchmark/LGASSERI_QUAKE_$timestamp.log 
echo "finished with code "$?

