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

mkdir -p ./ant_benchmark

echo "Starting E.coli is220 on ant08..."
rm -f ./ant_benchmark/ECOLI_IS220_QUAKE.log
srun -w ant08 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is220/ECOLI_IS220_QUAKE_BM.info  $1 >> ./ant_benchmark/ECOLI_IS220_QUAKE_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS220_QUAKE_$timestamp.log &
echo "done"

echo "Starting E.coli is480 on ant09..."
rm -f ./ant_benchmark/ECOLI_IS480_QUAKE.log
srun -w ant09 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is480/ECOLI_IS480_QUAKE_BM.info $1 >> ./ant_benchmark/ECOLI_IS480_QUAKE_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS480_QUAKE_$timestamp.log &
echo "done"

echo "Starting B.faecium on ant10..."
rm -f ./ant_benchmark/BFAECIUM_QAUKE.log
srun -w ant10 ./src/test/teamcity/teamcity.py /tmp/data/input/B.faecium/BFAECIUM_QAUKE_BM.info $1 >> ./ant_benchmark/BFAECIUM_QAUKE_$timestamp.log 2>> ./ant_benchmark/BFAECIUM_QAUKE_$timestamp.log &
echo "done"

echo "Starting L.gasseri on ant11..."
rm -f ./ant_benchmark/LGASSERI_QUAKE.log
srun -w ant11 ./src/test/teamcity/teamcity.py /tmp/data/input/L.gasseri/LGASSERI_QUAKE_BM.info $1 >> ./ant_benchmark/LGASSERI_QUAKE_$timestamp.log 2>> ./ant_benchmark/LGASSERI_QUAKE_$timestamp.log &
echo "done"

