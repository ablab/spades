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

mkdir -p ./ant_benchmark

echo "Starting E.coli is220 on ant16..."
rm -f ./ant_benchmark/ECOLI_IS220.log
srun -w ant16 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/is220/ECOLI_IS200.info >> ./ant_benchmark/ECOLI_IS220_$timestamp.log 2>> ./ant_benchmark/ECOLI_IS220_$timestamp.log &
echo "done"

echo "Starting E.coli UCSD lane 1 on ant17..."
rm -f ./ant_benchmark/ECOLI_SC_LANE_1.log
srun -w ant17 ./src/test/teamcity/teamcity.py /tmp/data/input/E.coli/sc_lane_1/ECOLI_SC_LANE_1.info >> ./ant_benchmark/ECOLI_SC_LANE_1_$timestamp.log 2>> ./ant_benchmark/ECOLI_SC_LANE_1_$timestamp.log &
echo "done"

echo "Starting M.ruber JGI lane 9 on ant18..."
rm -f ./ant_benchmark/MRUBER_JGI_LANE_9.log
srun -w ant18 ./src/test/teamcity/teamcity.py /tmp/data/input/M.ruber/jgi_lane_9/MRUBER_JGI_LANE_9.info >> ./ant_benchmark/MRUBER_JGI_LANE_9_$timestamp.log 2>> ./ant_benchmark/MRUBER_JGI_LANE_9_$timestamp.log &
echo "done"


echo "Starting P.heparinus JGI lane 7 on ant16..."
rm -f ./ant_benchmark/PHEPARINUS_JGI_LANE_7.log
srun -w ant16 ./src/test/teamcity/teamcity.py /tmp/data/input/P.heparinus/jgi_lane_7/PHEPARINUS_JGI_LANE_7.info >> ./ant_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log 2>> ./ant_benchmark/PHEPARINUS_JGI_LANE_7_$timestamp.log &
echo "done"

