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
srun -w ant08 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is220/ECOLI_IS220_QUAKE_BM.info $1 >> ./ant_benchmark/ECOLI_IS220_QUAKE_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS220_QUAKE_$timestamp.log &
echo "done"

echo "Starting E.coli UCSD lane 1 on ant09..."
rm -f ./ant_benchmark/ECOLI_SC_LANE_1.log
srun -w ant09 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/sc_lane_1/corrected/ECOLI_SC_LANE_1_BH_BM.info $1 >> ./ant_benchmark/ECOLI_SC_LANE_1_$timestamp.log 2>> ./ant_benchmark/ECOLI_SC_LANE_1_$timestamp.log &
echo "done"

echo "Starting M.ruber JGI lane 9 on ant10..."
rm -f ./ant_benchmark/MRUBER_JGI_LANE_9.log
srun -w ant10 ./src/test/teamcity/teamcity.py /tmp/data/input/M.ruber/jgi_lane_9/MRUBER_JGI_LANE_9_BM.info $1 >> ./ant_benchmark/MRUBER_JGI_LANE_9_$timestamp.log 2>> ./ant_benchmark/MRUBER_JGI_LANE_9_$timestamp.log &
echo "done"


echo "Starting P.heparinus JGI lane 7 on ant11..."
rm -f ./ant_benchmark/PHEPARINUS_JGI_LANE_7.log
srun -w ant11 ./src/test/teamcity/teamcity.py /tmp/data/input/P.heparinus/jgi_lane_7/PHEPARINUS_JGI_LANE_7_BM.info $1 >> ./ant_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log 2>> ./ant_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log &
echo "done"

