#!/bin/bash

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

echo "Starting E.coli is220 on ant16..."
rm -f ./ant_benchmark/ECOLI_IS220.log
srun -w ant16 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is220/ECOLI_IS200.info >> ./ant_benchmark/ECOLI_IS220_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS220_$timestamp.log &
echo "done"

echo "Starting E.coli is480 on ant17..."
rm -f ./ant_benchmark/ECOLI_IS480_QUAKE.log
srun -w ant17 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is480/ECOLI_IS480_QUAKE.info >> ./ant_benchmark/ECOLI_IS480_QUAKE_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS480_QUAKE_$timestamp.log &
echo "done"

echo "Starting L.grasseri on ant18..."
rm -f ./ant_benchmark/LGRASSERI_QUAKE.log
srun -w ant18 ./src/test/teamcity/teamcity.py /tmp/data/input/L.grasseri/LGRASSERI_QUAKE.info >> ./ant_benchmark/LGRASSERI_QUAKE_$timestamp.log 2>> ./ant_benchmark/LGRASSERI_QUAKE_$timestamp.log &
echo "done"


echo "Starting B.faecium on ant18..."
rm -f ./ant_benchmark/BFAECIUM_QAUKE.log
srun -w ant18 ./src/test/teamcity/teamcity.py /tmp/data/input/B.faecium/BFAECIUM_QAUKE.info >> ./ant_benchmark/BFAECIUM_QAUKE_$timestamp.log 2>> ./ant_benchmark/BFAECIUM_QAUKE_$timestamp.log &
echo "done"

