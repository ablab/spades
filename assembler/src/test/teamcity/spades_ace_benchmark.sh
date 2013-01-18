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

mkdir -p ./ace_benchmark

echo "Starting E.coli is220 on ace..."
rm -f ./ace_benchmark/ECOLI_IS220.log
./src/test/teamcity/teamcity.py /storage/acestorage/data/input/E.coli/is220/ECOLI_IS220_QUAKE_BM.info $1 >> ./ace_benchmark/ECOLI_IS220_$timestamp.log 2>> ./ace_benchmark/ECOLI_IS220_$timestamp.log
echo "finished with code "$?

echo "Starting E.coli UCSD lane 1 on ace..."
rm -f ./ace_benchmark/ECOLI_SC_LANE_1.log
./src/test/teamcity/teamcity.py /storage/acestorage/data/input/E.coli/sc_lane_1/bh20130116/ECOLI_SC_LANE_1_CORR.info $1 >> ./ace_benchmark/ECOLI_SC_LANE_1_$timestamp.log 2>> ./ace_benchmark/ECOLI_SC_LANE_1_$timestamp.log
echo "finished with code "$?

echo "Starting M.ruber JGI lane 9 on ace..."
rm -f ./ace_benchmark/MRUBER_JGI_LANE_9.log
./src/test/teamcity/teamcity.py /storage/acestorage/data/input/M.ruber/jgi_lane_9/MRUBER_JGI_LANE_9_MANUAL.info $1 >> ./ace_benchmark/MRUBER_JGI_LANE_9_$timestamp.log 2>> ./ace_benchmark/MRUBER_JGI_LANE_9_$timestamp.log
echo "finished with code "$?

echo "Starting P.heparinus JGI lane 7 on ace..."
rm -f ./ace_benchmark/PHEPARINUS_JGI_LANE_7.log
./src/test/teamcity/teamcity.py /storage/acestorage/data/input/P.heparinus/jgi_lane_7/PHEPARINUS_JGI_LANE_7_MANUAL.info $1 >> ./ace_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log 2>> ./ace_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log
echo "finished with code "$?

